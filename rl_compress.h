/* Original version of this file can be found here:
 * http://code.google.com/p/dcs-bwt-compressor/ */

/**********************************************************************
 * Byte stream compressor designed for compressing                    *
 * Burrows-Wheeler transforms.                                        *
 *                                                                    *
 * The main compression steps are the following:                      *
 * 1. Run length encoding: Replace a run of identical bytes bbb...bb  *
 *    with the pair (b, run_length).                                  *
 * 2. Turn each pair into a bit sequence: 8 bits of the byte followed *
 *    by the Elias gamma coding of the run length.                    *
 * 3. Encode each bit using range coding: Adaptive models are used    *
 *    for predicting the probabilities of bits.                       *
 **********************************************************************/

#ifndef DCSBWT_RL_COMPRESS_H__
#define DCSBWT_RL_COMPRESS_H__

#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include <cassert>

#include "stream_compressor.h"
#include "stream.h"
#include "globaldefs.h"

namespace dcsbwt {

extern int verbosity;

/* Probabilities are encoded as integers in [0,kProbabilityScale]
 * with p representing the probability p/kProbabilityScale. */
typedef int16 Probability;
static const int kLogProbabilityScale = 12;
static const Probability kProbabilityScale = (1 << kLogProbabilityScale);

/*********************************************************************
 * Entropy compressor for a sequence of bits.                        *
 * Given a probility distribution for each bit in the sequence,      *
 * ouputs a compressed representation of the sequence.               *
 * Uses range encoding (arithmetic encoding).                        *
 * The expected length of the compressed representation is           *
 * close to the information theoretic lower bound.                   *
 *********************************************************************/
class BitEncoder {
 public:
  BitEncoder() : low_(0), high_(0xFFFFFFFF) {}

  // Output is written to an OutStreamBuffer.
  void Connect(bwtc::OutStream* out) { output_ = out; }
  //TODO: Figure out what Disconnect should do if needed
  //bwtc::OutStream* Disconnect() { return output_.Disconnect(); }

  /* Append a bit to the sequence.
   * Even input cases bit==1(true), probability_of_one==0 and
   * bit==0(false), probability_of_one==kProbabilityScale are legal.
   * In these cases, up to four bytes of output is generated
   * from the single bit.  */
  void Encode(bool bit, Probability probability_of_one);

  /* Must be called to finish the encoding of a sequence.
   * After the call, BitEncoder is ready to start encoding a new sequence. */
  void Finish();

 private:
  uint32 low_;
  uint32 high_;
  bwtc::OutStream* output_;

  void EmitByte(unsigned char byte) { output_->WriteByte(byte); }
};

/*********************************************************************
 * Decompressor for a bit sequence compressed by BitEncoder.         *
 *********************************************************************/
class BitDecoder {
 public:
  BitDecoder() : low_(0), high_(0xFFFFFFFF) {}

  /* The compressed data is read from an InStreamBuffer. */
  void Connect(bwtc::InStream* in) { input_ = in; }
  //TODO: Do we need Disconnect()
  //bwtc::InStream* Disconnect() { return input_.Disconnect(); }

  /* Start() must be called to start the decoding of a sequence.
   * Nothing needs to be called to finish the decoding, but Start()
   * must be called again when starting to decode a new sequence.  */
  void Start();

  /* Get the next bit of the uncompressed sequence.
   * The probability distribution must be the same as the one used
   * when encoding the same bit. */
  bool Decode(Probability probability_of_one);

 private:
  uint32 low_;
  uint32 high_;
  uint32 next_;
  bwtc::InStream* input_;

  byte ReadByte() { return input_->ReadByte(); }
};



//////////////////////////////////////////////////////////////////
// BitPredictor maintains a probability distribution of a bit.
//
// The method Update(&probability, bit) adjusts the probability
// towards the given bit. This can be used for learning the probability.
//
// The amount of adjustment is controlled by two parameters:
// 1. update_delay: Higher value means a smaller adjustment.
//   REQUIRE: 0 <= update_delay <= kLogProbabilityScale
//   update_delay == 0 sets the probability to 0 or kProbabilityScale.
//   update_delay == kLogProbabilityScale means that the probability
//   changes by one (unless already at min or max).
// 2. min_probability limits the probability to the range
//   [min_probability, kProbabilityScale - min_probability]
//   Close to the extremes, adjustments are smaller.
//   REQUIRE: 0 <= min_probability < (2^update_delay)
//////////////////////////////////////////////////////////////////
class BitPredictor {
 public:
  BitPredictor(int update_delay, Probability min_probability) {
    SetParameters(update_delay, min_probability);
  }
  void SetParameters(int update_delay, Probability min_probability);
  void Update(Probability* probability_of_one, bool bit) const;
 private:
  int update_delay_;
  int rounding_correction_;
};

//////////////////////////////////////////////////////////////////
// BitHistory keeps an account of the distribution of
// the last few bits in a bit sequence.
//
// Bit history is essentially an automaton with a small number of states
// and two transitions (0-bit and 1-bit) out of each state.
// Half of the states represent situations where 0-bits are dominant
// in the recent history and the other half represent 1-dominant histories.
// There is no neutral state (except for history_length==0).
// Let's look at the 0-dominant states which we denote with
//  z1, z2, z3, ..., zk. 1-dominant states behave symmetrically.
// zi represents the bit sequence ...010101000...0,
// where the last run of zeros has length i.
// zk represents all such sequences with i >= k.
// The 0-transition from zi goes to z(i+1) (and from zk to itself).
// The 1-transition from zi goes either to zj for some j < i (for larger i)
// or to the first 1-dominant state (for smaller i).
// Here is an example:
//
// ----|    ---------------|
// |   v    |    |         v
// |--z4<--z3<--z2<--z1-->o1
//     |              ^
//     ---------------|
//
// The total number of states is 2^history_length.
// (REQUIRE: 0 <= history_length <= 5)
// The states are numbered from 0 to 2^history_length-1 in the
// order zk, ..., z1, o1, ..., ok.
// history_length==0 is a special case with only one state.
//
// A BitHistory object does not store a current state but supports
// computing the transitions using the method UpdateState().
// Thus the same BitHistory objet can be used with many bit sequences,
// each with its own current state.
//
// Each state can be associated with a probability of the next bit beeing 1.
// The probabilities can be updated according to observations using
// the method UpdateProbabilities(). This will update not only the probability
// of the given state but also those of its neighboring states
// (by a smaller amount). The updates are done using BitPredictors.
// Neighbor updates are not used if history_length <= 1.
//
// The state transitions and the probability updates are independent
// of each other, and there can be multiple independent probability
// vectors.
//////////////////////////////////////////////////////////////////
class BitHistory {
 public:
  typedef uint8 State;
  BitHistory(int history_length, int update_delay,
             Probability min_probability)
      : predictor_(update_delay, min_probability),
        neighbor_predictor_(update_delay + 1, min_probability) {
    SetHistoryLength(history_length);
  }

  void SetHistoryLength(int history_length);
  int HistoryLength() const { return history_length_; }
  int NumberOfStates() const { return 1 << history_length_; }
  State InitialState() const { return ((1 << history_length_) - 1) >> 1; }
  void UpdateState(State* state, bool bit) const;

  void SetPredictorParameters(int update_delay, Probability min_probability);
  // probability_vector must be a pointer to an array of size NumberOfStates()
  void UpdateProbability(Probability* probability_vector,
                         State state, bool bit) const;
  // adjustment changes initial probabilities from their normal values.
  // This is done by updating the probabilities abs(adjustment) times
  // towards zero if adjustment<0 or towards one if adjustment>0.
  void SetInitialProbabilities(Probability* probability_vector,
                               int adjustment = 0) const;
 private:
  State zero_transition_[32];
  State one_transition_[32];
  int history_length_;
  const Probability* initial_probabilities_;
  BitPredictor predictor_;
  BitPredictor neighbor_predictor_;
};

//////////////////////////////////////////////////////////////////
// BitSequencePredictorSet represents a set of adaptive predictors
// for bit sequences. The predictors are based on BitHistory.
// Each predictor is independent from others except that all are
// updated with the same BitHistory object, i.e.
// they have the same learning parameters.
//////////////////////////////////////////////////////////////////
class BitSequencePredictorSet {
 public:
  BitSequencePredictorSet(int num_sequences, int history_length,
                          int update_delay, Probability min_probability);

  // Returns the probability of the next bit being one.
  uint32 ProbabilityOfOne(int sequence) const;
  // Updates the prediction based on the actual next bit.
  void Update(int sequence, bool bit);

  static int64 SizeInBytes(int num_sequences, int history_length) {
    return 128 + num_sequences * sizeof(BitHistory::State)
        + (num_sequences << history_length) * sizeof(Probability);
  }
  int64 Size() const {
    return SizeInBytes(history_states_.size(), history_.HistoryLength());
  }
 private:
  int num_sequences_;
  int history_length_;
  BitHistory history_;
  std::vector<BitHistory::State> history_states_;
  std::vector<Probability> probabilities_;

  BitSequencePredictorSet(const BitSequencePredictorSet&);
  BitSequencePredictorSet& operator=(const BitSequencePredictorSet&);
};

//////////////////////////////////////////////////////////////////
// ByteCompressor is a simple byte stream compressor.
//
// Each byte is divided into bits and each bit has a context:
// the preceding bits in the byte. For example, if byte is
// 54=0011110, the 4th bit (a 1) has the context 001.
// Note that 001 is different from 1, 01, or 0001.
// There are 255 different contexts.
//
// The probability of each bit is predicted based on the previous
// bits in the same context. BitSequencePredictorSet is used
// for making the predictions for all contexts.
//
// The compressor has three parameters, which are the three
// parameters of BitHistory: history_length, update_delay,
// min_probability. The options string has one character,
// an ascii digit, for each parameter. The numeric value
// of the digit gives directly history_length, and update_delay,
// while the min_probability is computed using the formula
// min_probability = (2^digit)-1. For example,
// "344" means history_length=3, update_delay=4, min_probability=15.
//
// Legal parameter settings are all three digit combinations
// ("000".."999"), where the last digit is smaller or equal
// to the middle digit. For example "344" is legal but "345" is illegal.
// (The restriction comes from BitPredictor).
// An empty options string sets the parameters to default
// values (see ParseBitHistoryOptions in rl_compress.cc).
//////////////////////////////////////////////////////////////////
class ByteCompressor : public StreamCompressor {
 public:
  class Options : public StreamCompressor::OptionsBase {
   public:
    Options() { Set(""); }
    virtual bool Set(const std::string& options);
    virtual std::string Get() const;
    virtual int64 SizeInBytes() const {
      return 256 + BitSequencePredictorSet::SizeInBytes(256, parameters_[0]);
    }
    virtual StreamCompressor* GetCompressor() {
      return new ByteCompressor(parameters_);
    }
   private:
    int parameters_[3];
  };

  explicit ByteCompressor(const int* parameters)
      : predictor_(256, parameters[0], parameters[1],
                   (1 << parameters[2]) - 1) {
    if (verbosity > 2) {
      std::clog << "ByteCompressor"
                << " history_length=" << parameters[0]
                << " update_delay=" << parameters[1]
                << " min_probability=" << (1 << parameters[2]) - 1
                << std::endl;
    }
  }
  virtual ~ByteCompressor() {}

  virtual void ConnectPrivate(OutStreamBuffer* out) { encoder_.Connect(out); }
  virtual void DisconnectPrivate() { encoder_.Disconnect();}

  virtual void WriteBegin() {}
  virtual void WriteEndPrivate() { encoder_.Finish(); }
  virtual void Write(const char* bytes, size_t n) {
    for (; n; --n) { EncodeByte(*bytes++); }
  }

 private:
  BitSequencePredictorSet predictor_;
  BitEncoder encoder_;

  void EncodeByte(unsigned char byte);

  ByteCompressor(const ByteCompressor&);
  ByteCompressor& operator=(const ByteCompressor&);
};

//////////////////////////////////////////////////////////////////
// ByteDecompressor decompresses data compressed by ByteCompressor
//////////////////////////////////////////////////////////////////
class ByteDecompressor : public StreamDecompressor {
 public:
  static StreamDecompressor* Create(const std::string& options);
  explicit ByteDecompressor(const int* parameters)
      : predictor_(256, parameters[0], parameters[1],
                   (1 << parameters[2]) - 1) {
    if (verbosity > 2) {
      std::clog << "ByteDecompressor"
                << " history_length=" << parameters[0]
                << " update_delay=" << parameters[1]
                << " min_probability=" << (1 << parameters[2]) - 1
                << std::endl;
    }
  }
  virtual ~ByteDecompressor() {}

  virtual void ConnectPrivate(InStreamBuffer* out) { decoder_.Connect(out); }
  virtual void DisconnectPrivate() { decoder_.Disconnect(); }

  virtual void ReadBegin() { decoder_.Start(); }
  virtual void ReadEnd() {}
  virtual void Read(char* bytes, size_t n) {
    for (; n; --n) *bytes++ = DecodeByte();
  }

  virtual int64 SizeInBytes() const {
    return predictor_.Size();
  }

 private:
  BitSequencePredictorSet predictor_;
  BitDecoder decoder_;

  unsigned char DecodeByte();

  ByteDecompressor(const ByteDecompressor&);
  ByteDecompressor& operator=(const ByteDecompressor&);
};


//////////////////////////////////////////////////////////////////
// CharCoder compresses and decompresses byte streams.
//
// Differences to ByteCompressor/ByteDecompressor:
// - CharCoder is not a selfstanding StreamCompressor but
//   is used from within RunLengthCompressor and RunLengthDecompressor.
// - CharCoder is particularly designed for compressing
//   byte streams resulting from run length encoding.
// - CharCoder uses a more complex predictor model.
// - CharCoder does both compressing and decompressing
//   (in the same private function) to ensure consistency
//   under the complex predictor.
//
// Like ByteCompressor, the prediction model looks each bit in
// the context of preceding bits in the byte (the bit context)
// and uses a separate BitHistory state for each bit context.
// The differences are:
// - BitHistory state updates are delayed by one byte.
//   The rational is that after run length encoding
//   consecutive bytes are always different. Thus,
//   the preceding byte is used for updating the BitHistory states
//   only after the update can no more affect the encoding
//   of the current byte.
// - In addition to bit context, there is a byte context,
//   which depends on the preceding 0-3 _distinct_ bytes.
//   For each preceding bytes, the model takes into account,
//   whether it has the same bit context, and if it has
//   what the bit in the current bit position is.
//   Thus there are 1, 3, 9 or 27 different byte contexts.
//
// The prediction is derived as follows:
// - Use the combination of the bit context and the byte context
//   to choose a BitHistory probability vector P.
// - Use the bit context alone to choose a BitHistory state s.
// - The prediction P[s].
// Updating gos as follows:
// - P[s] (and its neighbors P[s-1] and P[s+1]) are updated
//   immediately after getting the prediction.
// - The BitHistory state is updated in the delayed manner
//   explained above.
//
// The parameters are given as an array of 4 integers. They are:
// 0: byte_context_length, the number of bytes used for byte context
// 1: history_length
// 2: update_delay
// 3: log_2(min_probability + 1)
// The last three are the same as the parameters of ByteCompressor
//
// TODO: Allow a varying length encoding of bytes, such as
// Huffmann encoding, instead of the even 8-bit encoding.
// The encoding must be static; the prediction mechanism assumes
// that the same byte is always encoded in the same way.
// The main purpose of this change would be to speed up coding
// by reducing the number bits to code. It might improve compression, too.
//////////////////////////////////////////////////////////////////
class CharCoder {
 public:
  explicit CharCoder(const int* parameters);
  void SetEncoder(BitEncoder* encoder) {
    assert(decoder_ == NULL);
    encoder_ = encoder;
  }
  void SetDecoder(BitDecoder* decoder) {
    assert(encoder_ == NULL);
    decoder_ = decoder;
  }
  void Encode(unsigned char byte) {
    assert(encoder_ != NULL);
    Code(&byte);
  }
  unsigned char Decode() {
    assert(decoder_ != NULL);
    unsigned char byte = 0;
    Code(&byte);
    return byte;
  }

  static int64 SizeInBytes(const int* parameters) {
    int byte_context_length = parameters[0];
    int history_length = parameters[1];
    int num_probabilities = 256 << (history_length + 2 * byte_context_length);
    return 128 + 256 * sizeof(BitHistory::State)
        + num_probabilities * sizeof(Probability);
  }
  int64 Size() const {
    int parameters[2];
    parameters[0] = byte_context_length_;
    parameters[1] = history_.HistoryLength();
    return SizeInBytes(parameters);
  }
 private:
  int byte_context_length_;
  unsigned char recent_bytes_[3];
  BitHistory history_;
  std::vector<BitHistory::State> history_states_;
  std::vector<Probability> probabilities_;
  BitEncoder* encoder_;
  BitDecoder* decoder_;

  void Code(unsigned char* byte);
};

//////////////////////////////////////////////////////////////////
// RunLengthCoder compresses and decompresses run lengths in
// run length encoding.
//
// The encoding has two steps:
// 1. Turn the run length into a bit sequence by Elias gamma code.
// 2. Entropy code each bit using a prediction model based on BitHistory
//    (and using BitSequencePredictorSet to represent the models).
//
// Elias gamma code has two parts:
// 1. Unary encoding of the number of significant bits (to be precise,
//    the number of bits after the most significant 1-bit).
// 2. The significant bits (without the most significant 1-bit).
// For example:
//  1 =    ...1 ->    0 (no bits after the most sigificant 1-bit)
//  2 =   ...10 ->   10 0
//  3 =   ...11 ->   10 1
//  4 =  ...100 ->  110 00
// ...
// 13 = ...1101 -> 1110 101
// (Side note: Common description writes the unary part with
//  the roles of 0s and 1s reversed. Then the last 1-bit of the unary
//  part can also be seen as the omitted most significant 1-bit.
//  On the other hand, under the present encoding, the bit sequences
//  have the correct lexicographic order.)
//
// There are three kinds of prediction models:
// 1. The first bit (which tells whether the run length is one or more)
//    is predicted in the context of the associated byte.
//    I.e., there is a separate predictor for each byte.
// 2. The rest of the bytes are predicted without context.
//    I.e., there is a separate predictor for each position in the unary part.
// 3. The significant bits are predicted in the context of
//    the number of significant bits. I.e., there is a predictor
//    for each combination of the number of significant bits and
//    the bit positions in the significant bits.
//
// Each of the three prediction model kinds has its own
// BitSequencePredictorSet object and thus can have its
// own learning parameters. These are provided in an array
// of 9 integers, three for each kind (in the above order).
// The three parameters for one kind are:
// 1: history_length
// 2: update_delay
// 3: log_2(min_probability + 1)
//////////////////////////////////////////////////////////////////
class RunLengthCoder {
 public:
  RunLengthCoder(int log_max_run_length, const int* parameters);
  void SetEncoder(BitEncoder* encoder) {
    assert(decoder_ == NULL);
    encoder_ = encoder;
  }
  void SetDecoder(BitDecoder* decoder) {
    assert(encoder_ == NULL);
    decoder_ = decoder;
  }

  void Encode(unsigned char byte, int64 run_length) {
    assert(encoder_ != NULL);
    assert(run_length >= 1);
    Code(byte, &run_length);
  }
  int64 Decode(unsigned char byte) {
    assert(decoder_ != NULL);
    int64 run_length = 0;
    Code(byte, &run_length);
    assert(run_length >= 1);
    return run_length;
  }

  static int64 SizeInBytes(int log_max_run_length, const int* parameters) {
    int num_contexts = (log_max_run_length * (log_max_run_length + 1)) >> 1;
    return 32
        + BitSequencePredictorSet::SizeInBytes(256, parameters[0])
        + BitSequencePredictorSet::SizeInBytes(log_max_run_length,
                                               parameters[3])
        + BitSequencePredictorSet::SizeInBytes(num_contexts,
                                               parameters[6]);
  }
  int64 Size() const {
    return 32 + first_bit_predictor_.Size()
        + unary_part_predictor_.Size() + significant_bits_predictor_.Size();
  }
 private:
  int64 max_run_length_;
  BitSequencePredictorSet first_bit_predictor_;
  BitSequencePredictorSet unary_part_predictor_;
  BitSequencePredictorSet significant_bits_predictor_;
  BitEncoder* encoder_;
  BitDecoder* decoder_;

  void Code(unsigned char byte, int64* run_length);
};

//////////////////////////////////////////////////////////////////
// RunLengthCompressor is a byte stream compressor
// based on run length encoding followed by entropy encoding.
//
// The character of each run is encoded with CharCoder
// and the run length with RunLengthCoder.
//////////////////////////////////////////////////////////////////
class RunLengthCompressor : public StreamCompressor {
 public:

  class Options : public StreamCompressor::OptionsBase {
   public:
    static const int kNumParameters = 13;
    static const int kLogMaxRunLength = 31; // max run length = 2^32 - 1
    Options() { Set(""); }
    virtual bool Set(const std::string& options);
    virtual std::string Get() const;
    virtual int64 SizeInBytes() const {
      return 1024 + CharCoder::SizeInBytes(parameters_)
          + RunLengthCoder::SizeInBytes(kLogMaxRunLength, parameters_ + 4)
          + 4800;  // statistics
    }
    virtual StreamCompressor* GetCompressor() {
      if (verbosity > 1) {
        std::clog << "Using RunLengthCompressor with options=" << Get()
                  << std::endl;
      }
      return new RunLengthCompressor(parameters_);
    }
   private:
    int parameters_[kNumParameters];
  };

  explicit RunLengthCompressor(const int *parameters)
      : char_coder_(parameters),
        run_length_coder_(Options::kLogMaxRunLength, parameters + 4) {}
  virtual ~RunLengthCompressor() {}

  virtual void ConnectPrivate(OutStreamBuffer* out) {
    statistics_collector.StartEncoding();
    encoder_.Connect(out);
    char_coder_.SetEncoder(&encoder_);
    run_length_coder_.SetEncoder(&encoder_);
  }
  virtual void DisconnectPrivate() {
    run_length_coder_.SetEncoder(NULL);
    char_coder_.SetEncoder(NULL);
    encoder_.Disconnect();
    statistics_collector.EndEncoding();
  }

  virtual void WriteBegin() {
    current_run_length_ = 0;
    statistics_collector.StartBlock();
  }
  virtual void WriteEndPrivate() {
    EncodeRun(current_char_, current_run_length_);
    current_run_length_ = 0;
    encoder_.Finish();
    statistics_collector.EndBlock();
  }
  virtual void Write(const char* bytes, size_t n);

 private:
  CharCoder char_coder_;
  RunLengthCoder run_length_coder_;
  BitEncoder encoder_;
  unsigned char current_char_;
  int64 current_run_length_;

  void EncodeRun(unsigned char byte, int64 run_length);

  class StatisticsCollector {
   public:
    void StartEncoding();
    void EndEncoding();
    void StartBlock();
    void EndBlock();
    void RecordRun(unsigned char byte, int64 run_length);
   private:
    static const int kLogLongRunLength = 6;
    static const int kLongRunLength = 1 << kLogLongRunLength;
    int64 num_chars_[256];
    int64 num_runs_by_char_[256];
    int64 num_short_runs_by_length_[kLongRunLength];
    int64 num_long_runs_by_log_length_[64 - kLogLongRunLength];
  };
  StatisticsCollector statistics_collector;

  RunLengthCompressor(const RunLengthCompressor&);
  RunLengthCompressor& operator=(const RunLengthCompressor&);
};

//////////////////////////////////////////////////////////////////
// RunLengthDecompressor decompresses data compressed with
// RunLengthCompressor
//////////////////////////////////////////////////////////////////
class RunLengthDecompressor : public StreamDecompressor {
 public:
  static const int kNumParameters = 13;
  static const int kLogMaxRunLength = 31; // max run length = 2^32 - 1
  static StreamDecompressor* Create(const std::string& options);
  explicit RunLengthDecompressor(const int* parameters)
      : char_coder_(parameters),
        run_length_coder_(kLogMaxRunLength, parameters + 4) {}
  virtual ~RunLengthDecompressor() {}

  virtual void ConnectPrivate(InStreamBuffer* out) {
    decoder_.Connect(out);
    char_coder_.SetDecoder(&decoder_);
    run_length_coder_.SetDecoder(&decoder_);
  }
  virtual void DisconnectPrivate() {
    run_length_coder_.SetDecoder(NULL);
    char_coder_.SetDecoder(NULL);
    decoder_.Disconnect();
  }

  virtual void ReadBegin() {
    decoder_.Start();
    current_run_length_ = 0;
  }
  virtual void ReadEnd() { current_run_length_ = 0; }
  virtual void Read(char* bytes, size_t n);

  virtual int64 SizeInBytes() const {
    return 1024 + char_coder_.Size() + run_length_coder_.Size();
  }

 private:
  CharCoder char_coder_;
  RunLengthCoder run_length_coder_;
  BitDecoder decoder_;
  unsigned char current_char_;
  int64 current_run_length_;

  void DecodeRun(unsigned char* byte, int64* run_length);

  RunLengthDecompressor(const RunLengthDecompressor&);
  RunLengthDecompressor& operator=(const RunLengthDecompressor&);
};

}  // namespace dcsbwt

#endif  // DCSBWT_RL_COMPRESS_H__
