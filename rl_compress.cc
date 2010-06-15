/* Original implementations for BitEncoder and Decoder can be found here: *
 * http://code.google.com/p/dcs-bwt-compressor/                           */

#include <cassert>

#include <string>

#include "rl_compress.h"
#include "stream.h"
#include "globaldefs.h" /* Important definitions */

namespace dcsbwt {

extern int statistics;

/**************************************************************************
 * BitEncoder and BitDecoder maintain a range [low_,high_] (of uint32's). *
 * The next four bytes of the compressed output/input represent           *
 * some value in this range.                                              *
 * For each bit (en/de)coded, the range is split into two parts           *
 * with the sizes proportional to the probabilities of 1-bit and 0-bit.   *
 * BitEncoder chooses one the two subranges as the new range [low_,high_].*
 * based on the actual value of the bit.                                  *
 * BitDecoder determines the bit by finding out, which subrange           *
 * contains the value represented by the next four input bytes.           *
 * Whenever low_ and high_ share the most significant byte,               *
 * the range is expanded by the factor 256 shifting the shared byte out.  *
 * BitEncoder emits the shared byte. BitDecoder shifts the shared byte    *
 * out of the next four input bytes that it is holding and reads in       *
 * a new input byte.                                                      *
 **************************************************************************/

namespace {
/* Split a range [low,high] into [low,split] and [split+1,high]
 * proportional to the probability and its complement. */
uint32 Split(uint32 low, uint32 high, uint32 probability) {
  assert(probability <= kProbabilityScale);
  assert(low < high);
  /* range_size is high-low-1 rather than high-low+1 to ensure that
   * neither subrange is empty, even when probability is 0 and 1. */
  uint32 range_size = high - low - 1;
  /* split = low + round(range_size * p)
   *       = low + floor(range_size * p + .5)
   * where p is the probability as a real value. */
  uint32 high_bits = range_size >> kLogProbabilityScale;
  uint32 low_bits = range_size & (kProbabilityScale - 1);
  uint32 half = kProbabilityScale >> 1;
  uint32 split = low + high_bits * probability
      + ((low_bits * probability + half) >> kLogProbabilityScale);
  assert(split >= low);
  assert(split < high);
  return split;
}
}  // namespace

BitEncoder::BitEncoder()
    : low_(0), high_(0xFFFFFFFF), counter_(0), output_(NULL) {}

BitEncoder::~BitEncoder() { }

void BitEncoder::Encode(bool bit, Probability probability_of_one) {
  if (verbosity > 7) {
    std::clog << "Encoding bit " << int(bit) << " with probability "
              << (bit ? probability_of_one :
                  kProbabilityScale - probability_of_one)
              << std::endl;
  }
  uint32 split = Split(low_, high_, probability_of_one);
  /* Choose the subrange. */
  if (bit) high_ = split; else low_ = split + 1;
  while (((low_ ^ high_) & 0xFF000000) == 0) {
    EmitByte(low_ >> 24);
    low_ <<= 8;
    high_ = (high_ << 8) + 255;
  }
  assert(low_ < high_);
}

void BitEncoder::Finish() {
  /* Emit 4 bytes representing any value in [low_,high_]
   * These are needed for BitDecoder lookahead. */
  EmitByte(low_ >> 24);
  EmitByte(255);
  EmitByte(255);
  EmitByte(255);
  output_->Flush();
  /* Prepare to encode another sequence. */
  low_ = 0;
  high_ = 0xFFFFFFFF;
}

BitDecoder::BitDecoder() :
    low_(0), high_(0xFFFFFFFF), next_(0), input_(NULL) {}

BitDecoder::~BitDecoder() { }

void BitDecoder::Start() {
  low_ = 0;
  high_ = 0xFFFFFFFF;
  next_ = ReadByte();
  next_ = (next_ << 8) + ReadByte();
  next_ = (next_ << 8) + ReadByte();
  next_ = (next_ << 8) + ReadByte();
}

bool BitDecoder::Decode(Probability probability_of_one) {
  uint32 split = Split(low_, high_, probability_of_one);
  bool bit = (next_ <= split);
  if (bit) high_ = split; else low_ = split + 1;
  while (((low_ ^ high_) & 0xFF000000) == 0) {
    low_ <<= 8;
    high_ = (high_ << 8) + 255;
    next_ = (next_ << 8) + ReadByte();
  }
  assert(low_ < high_);
  assert(next_ >= low_);
  assert(next_ <= high_);
  if (verbosity > 7) {
    std::clog << "Decoded bit " << int(bit) << " with probability "
              << (bit ? probability_of_one :
                  kProbabilityScale - probability_of_one)
              << std::endl;
  }
  return bit;
}

#if 0

//////////////////////////////////////////////////////////////////
// BitPredictor
//
// Note: The restriction min_probability < 2^update_delay
// comes from the requirement that the expressions (in Update)
//   kProbabilityScale - *probability + rounding_correction_
//   *probability + rounding_correction_
// need to be non-negative. This is because the result of
// the shift (>> update_delay_) for negative values
// is implementation-defined (5.8/3).
//////////////////////////////////////////////////////////////////
void BitPredictor::SetParameters(int update_delay,
                                 Probability min_probability) {
  assert(update_delay >= 0);
  assert(update_delay <= kLogProbabilityScale);
  assert(min_probability >= 0);
  assert(min_probability < (1 << update_delay));
  update_delay_ = update_delay;
  rounding_correction_ = (1 << update_delay) - 1 - min_probability;
  assert(rounding_correction_ >= 0);
  assert(rounding_correction_ < 1 << update_delay);
}
void BitPredictor::Update(Probability* probability, bool bit) const {
  assert(*probability >= 0);
  assert(*probability <= kProbabilityScale);
  if (bit)
    *probability += ((kProbabilityScale - *probability + rounding_correction_)
                     >> update_delay_);
  else
    *probability -= ((*probability + rounding_correction_) >> update_delay_);
  assert(*probability >= 0);
  assert(*probability <= kProbabilityScale);
}


//////////////////////////////////////////////////////////////////
// BitHistory
//////////////////////////////////////////////////////////////////
namespace {
// BitHistory state transitions.
// Only 1-transitions of 0-dominant states are given.
// 0-transitions of 0-dominant states always go to the previous state.
// The transitions of 1-dominant states are symmetric.
const BitHistory::State bit_history_transition_1[1] = {1};
const BitHistory::State bit_history_transition_2[2] = {2, 2};
const BitHistory::State bit_history_transition_3[4] = {3, 4, 4, 4};
const BitHistory::State bit_history_transition_4[8] = {6, 6, 6, 7, 7, 7, 8, 8};
const BitHistory::State bit_history_transition_5[16] = {
  11, 12, 12, 12, 12, 12, 12, 13, 13, 13, 14, 14, 15, 15, 16, 16
};
const BitHistory::State* bit_history_transitions[6] = {
  NULL,
  bit_history_transition_1,
  bit_history_transition_2,
  bit_history_transition_3,
  bit_history_transition_4,
  bit_history_transition_5
};

// BitHistory initial probabilities.
const Probability bit_history_probability_0[1] = {2048};
const Probability bit_history_probability_1[1] = {1358};
const Probability bit_history_probability_2[2] = {768, 1408};
const Probability bit_history_probability_3[4] = {636, 864, 1196, 1675};
const Probability bit_history_probability_4[8] = {
  391, 469, 572, 704, 876, 1100, 1391, 1768
};
const Probability bit_history_probability_5[16] = {
  214, 233, 256, 284, 318, 360, 411, 474,
  550, 642, 755, 893, 1061, 1266, 1516, 1821
};
const Probability* bit_history_probabilities[6] = {
  bit_history_probability_0,
  bit_history_probability_1,
  bit_history_probability_2,
  bit_history_probability_3,
  bit_history_probability_4,
  bit_history_probability_5
};
}

void BitHistory::SetHistoryLength(int history_length) {
  assert(history_length >= 0);
  assert(history_length <= 5);
  history_length_ = history_length;
  int num_states = (1 << history_length);
  initial_probabilities_ = bit_history_probabilities[history_length];
  const State* transitions = bit_history_transitions[history_length];
  for (int i = 0; i < (num_states >> 1); ++i) {
    zero_transition_[i] = i - 1;
    one_transition_[i] = transitions[i];
    zero_transition_[num_states - i - 1] = num_states - transitions[i] - 1;
    one_transition_[num_states - i - 1] = num_states - i;
  }
  zero_transition_[0] = 0;
  one_transition_[num_states - 1] = num_states - 1;
  for (int i = 0; i < num_states; ++i) {
    assert(int(zero_transition_[i]) < num_states);
    assert(int(one_transition_[i]) < num_states);
  }
}
void BitHistory::UpdateState(State* state, bool bit) const {
  assert(int(*state) < (1 << history_length_));
  if (bit) *state = one_transition_[*state];
  else *state = zero_transition_[*state];
  assert(int(*state) < (1 << history_length_));
}

void BitHistory::SetPredictorParameters(int update_delay,
                                        Probability min_probability) {
  predictor_.SetParameters(update_delay, min_probability);
  int neighbor_update_delay = std::min(kLogProbabilityScale, update_delay + 1);
  neighbor_predictor_.SetParameters(neighbor_update_delay, min_probability);
}
void BitHistory::SetInitialProbabilities(Probability* probability_vector,
                                         int adjustment) const {
  int num_states = (1 << history_length_);
  for (int state = 0; state < num_states; ++ state) {
    if (state < (num_states >> 1)) {
      probability_vector[state] = initial_probabilities_[state];
    } else {
      probability_vector[state]
          = initial_probabilities_[num_states - state - 1];
    }
    while (adjustment > 0) {
      predictor_.Update(&probability_vector[state], true);
      --adjustment;
    }
    while (adjustment < 0) {
      predictor_.Update(&probability_vector[state], false);
      ++adjustment;
    }
  }
}
void BitHistory::UpdateProbability(Probability* probability_vector,
                         State state, bool bit) const {
  int num_states = (1 << history_length_);
  //assert(int(state) >= 0);
  assert(int(state) < num_states);
  predictor_.Update(&probability_vector[state], bit);
  if (history_length_ > 2) {
    if (state > 0)
      neighbor_predictor_.Update(&probability_vector[state - 1], bit);
    if (state < num_states - 1)
      neighbor_predictor_.Update(&probability_vector[state + 1], bit);
  }
}

//////////////////////////////////////////////////////////////////
// BitSequencePredictorSet
//////////////////////////////////////////////////////////////////
BitSequencePredictorSet::BitSequencePredictorSet(
    int num_sequences, int history_length,
    int update_delay, Probability min_probability)
    : num_sequences_(num_sequences),
      history_(history_length, update_delay, min_probability),
      history_states_(num_sequences),
      probabilities_(num_sequences << history_length) {
  for (int sequence = 0; sequence < num_sequences; ++sequence) {
    history_states_[sequence] = history_.InitialState();
    history_.SetInitialProbabilities(
        &probabilities_[sequence << history_length]);
  }
}
uint32 BitSequencePredictorSet::ProbabilityOfOne(int sequence) const {
  assert(sequence >= 0);
  assert(sequence <= num_sequences_);
  BitHistory::State state = history_states_[sequence];
  int history_length = history_.HistoryLength();
  int position = (sequence << history_length) + state;
  assert(position >= 0);
  assert(position < probabilities_.size());
  return probabilities_[position];
}
void BitSequencePredictorSet::Update(int sequence, bool bit) {
  assert(sequence >= 0);
  assert(sequence < num_sequences_);
  BitHistory::State* state = &history_states_[sequence];
  int history_length = history_.HistoryLength();
  history_.UpdateProbability(&probabilities_[sequence << history_length],
                             *state, bit);
  history_.UpdateState(state, bit);
}


//////////////////////////////////////////////////////////////////
// ByteCompressor and ByteDecompressor
//////////////////////////////////////////////////////////////////
namespace {
// Each context is turned into a unique number by prepending a 1-bit.
// For example 001 -> 1001=9 and 01 -> 101=5
// These functions extract the bit and the context
// given a byte and a bit position.
// They are used outside ByteCompressor and ByteDecompressor, too,
// in CharCoder, for example.
int GetBit(int byte, int bit_position) {
  return (byte >> bit_position) & 1;
}
int BitContext(int byte, int bit_position) {
  return (byte | 256) >> (bit_position + 1);
}

// Transform a string representation of BitHistory options
// into an array of integers.
bool ParseBitHistoryOptions(const std::string& options,
                            int* parameters) {
  if (options.length() != 3) return false;
  int history_length = options[0] - '0';
  int update_delay = options[1] - '0';
  int log_min_probability = options[2] - '0';
  if (verbosity > 4) {
    std::clog << "Parsing BitHistory options"
              << " history_length=" << history_length
              << " update_delay=" << update_delay
              << " log_min_probability=" << log_min_probability
              << std::endl;
  }
  if (history_length < 0 || history_length > 9) return false;
  if (update_delay < 0 || update_delay > 9) return false;
  if (log_min_probability < 0 || log_min_probability > update_delay)
    return false;
  parameters[0] = history_length;
  parameters[1] = update_delay;
  parameters[2] = log_min_probability;
  return true;
}
}

bool ByteCompressor::Options::Set(const std::string& options) {
  std::string opts = options;
  if (opts.length() == 0) opts = "344";
  return ParseBitHistoryOptions(opts, parameters_);
}
std::string ByteCompressor::Options::Get() const {
  std::string options;
  for (int i = 0; i < 3; ++i) {
    assert(parameters_[i] >= 0);
    assert(parameters_[i] <= 9);
    options += static_cast<char>('0' + parameters_[i]);
  }
  return options;
}
void ByteCompressor::EncodeByte(unsigned char byte) {
  if (verbosity > 6) {
    std::clog << "ByteCompressor: encoding byte " << int(byte) << std::endl;
  }
  for (int i = 7; i >= 0; --i) {
    int bit = GetBit(byte, i);
    int context = BitContext(byte, i);
    encoder_.Encode(bit, predictor_.ProbabilityOfOne(context));
    predictor_.Update(context, bit);
  }
}

StreamDecompressor* ByteDecompressor::Create(const std::string& options) {
  int parameters[3];
  bool ok = ParseBitHistoryOptions(options, parameters);
  if (!ok) return NULL;
  return new ByteDecompressor(parameters);
}
unsigned char ByteDecompressor::DecodeByte() {
  int context = 1;
  for (int i = 7; i >= 0; --i) {
    bool bit = decoder_.Decode(predictor_.ProbabilityOfOne(context));
    predictor_.Update(context, bit);
    context = (context << 1) + bit;
  }
  unsigned char byte = context & 255;
  if (verbosity > 6) {
    std::clog << "ByteDecompressor: decoded byte " << int(byte) << std::endl;
  }
  return byte;
}

//////////////////////////////////////////////////////////////////
// CharCoder
//////////////////////////////////////////////////////////////////
namespace {
// Get a bit if in given context. The return value is one of:
//  0: same context, bit is 0
//  1: same context, bit is 1
//  2: different context
static int ByteContext(int byte, int bit_position, int other_context) {
  int context = BitContext(byte, bit_position);
  return context == other_context ? GetBit(byte, bit_position) : 2;
}
}

CharCoder::CharCoder(const int* parameters)
    : byte_context_length_(parameters[0]),
      history_(parameters[1], parameters[2],
               (1 << parameters[3]) - 1),
      history_states_(256),
      probabilities_(256 << (history_.HistoryLength()
                             + (2 * byte_context_length_))),
      encoder_(NULL), decoder_(NULL) {
  if (verbosity > 3) {
    std::clog << "CharCoder byte_context_length=" << byte_context_length_
              << " history_length=" << history_.HistoryLength()
              << " update_delay=" << parameters[2]
              << " min_probability=" << ((1 << parameters[3]) - 1)
              << std::endl;
  }
  for (int i = 0; i < 3; ++i) { recent_bytes_[i] = '\0'; }
  for (int bit_context = 0; bit_context < 256; ++bit_context) {
    history_states_[bit_context] = history_.InitialState();
    int byte_context_begin = bit_context << (2 * byte_context_length_);
    int byte_context_end = (bit_context + 1) << (2 * byte_context_length_);
    for (int byte_context = byte_context_begin;
         byte_context < byte_context_end; ++byte_context) {
      history_.SetInitialProbabilities(
          &probabilities_[byte_context << history_.HistoryLength()]);
    }
  }
}

void CharCoder::Code(unsigned char* byte) {
  int history_length = history_.HistoryLength();
  int bit_context = 1;
  for (int i = 7; i >= 0; --i) {
    int byte_context = 0;
    for (int j = 0; j < byte_context_length_; ++j) {
      byte_context <<= 2;
      byte_context += ByteContext(recent_bytes_[j], i, bit_context);
    }
    int full_context
        = (bit_context << (byte_context_length_ << 1)) + byte_context;
    BitHistory::State state = history_states_[bit_context];
    int position = (full_context << history_length) + state;
    Probability probability = probabilities_[position];
    bool bit;
    if (encoder_ != NULL) {
      bit = GetBit(*byte, i);
      encoder_->Encode(bit, probability);
    } else {
      bit = decoder_->Decode(probability);
    }
    history_.UpdateProbability(
        &probabilities_[full_context << history_length], state, bit);
    history_.UpdateState(&history_states_[BitContext(recent_bytes_[0], i)],
                         GetBit(recent_bytes_[0], i));
    bit_context = (bit_context << 1) + bit;
  }
  *byte = (bit_context & 255);
  // Update the last three distinct bytes.
  // This relies on the post-run-length-encoding property
  // that consecutive bytes are always different.
  // Thus we always have the last two bytes but the third byte
  // might be more distant.
  if (*byte != recent_bytes_[1]) {
    recent_bytes_[2] = recent_bytes_[1];
  }
  recent_bytes_[1] = recent_bytes_[0];
  recent_bytes_[0] = *byte;
}

//////////////////////////////////////////////////////////////////
// RunLengthCoder
//////////////////////////////////////////////////////////////////
namespace {
// Return the number of context when the nunber of bits
// is smaller or equal to num_bits.
// This is simply 1+2+3+...+num_bits
int SignificantBitContexts(int num_bits) {
  return (num_bits * (num_bits + 1)) >> 1;
}
}
RunLengthCoder::RunLengthCoder(int log_max_run_length, const int* parameters)
    : max_run_length_((1LL << (log_max_run_length + 1)) - 1),
      first_bit_predictor_(256, parameters[0], parameters[1],
                           (1 << parameters[2]) - 1),
      unary_part_predictor_(log_max_run_length, parameters[3], parameters[4],
                            (1 << parameters[5]) - 1),
      significant_bits_predictor_(
          SignificantBitContexts(log_max_run_length),
          parameters[6], parameters[7], (1 << parameters[8]) - 1),
      encoder_(NULL), decoder_(NULL) {
  assert(log_max_run_length <= 62);
  if (verbosity > 3) {
    std::clog << "RunlengthCoder max_run_length=" << max_run_length_
              << std::endl;
    std::clog << "RunlengthCoder first_bit_predictor:"
              << " history_length=" << parameters[0]
              << " update_delay=" << parameters[1]
              << " min_probability=" << ((1 << parameters[2]) - 1)
              << std::endl;
    std::clog << "RunlengthCoder unary_part_predictor:"
              << " history_length=" << parameters[3]
              << " update_delay=" << parameters[4]
              << " min_probability=" << ((1 << parameters[5]) - 1)
              << std::endl;
    std::clog << "RunlengthCoder significant_bits_predictor:"
              << " history_length=" << parameters[6]
              << " update_delay=" << parameters[7]
              << " min_probability=" << ((1 << parameters[8]) - 1)
              << std::endl;
  }
}

void RunLengthCoder::Code(unsigned char byte, int64* run_length) {
  assert(*run_length <= max_run_length_);
  bool bit;

  // Code the first bit
  Probability probability = first_bit_predictor_.ProbabilityOfOne(byte);
  if (encoder_ != NULL) {
    bit = (1 < *run_length);
    encoder_->Encode(bit, probability);
  } else {
    assert(decoder_ != NULL);
    bit = decoder_->Decode(probability);
  }
  first_bit_predictor_.Update(byte, bit);

  // Code the rest of the unary part
  int bits_after_leading_one = 0;
  while(bit) {
    ++bits_after_leading_one;
    probability
        = unary_part_predictor_.ProbabilityOfOne(bits_after_leading_one);
    if (encoder_ != NULL) {
      bit = ((*run_length >> bits_after_leading_one) > 1);
      encoder_->Encode(bit, probability);
    } else {
      assert(decoder_ != NULL);
      bit = decoder_->Decode(probability);
    }
    unary_part_predictor_.Update(bits_after_leading_one, bit);
  }

  // Code the significant bits
  int context_begin = SignificantBitContexts(bits_after_leading_one - 1);
  int64 length = 1;
  for (int i = bits_after_leading_one - 1; i >= 0; --i) {
    int context = context_begin + i;
    probability = significant_bits_predictor_.ProbabilityOfOne(context);
    if (encoder_ != NULL) {
      bit = (*run_length >> i) & 1;
      encoder_->Encode(bit, probability);
    } else {
      assert(decoder_ != NULL);
      bit = decoder_->Decode(probability);
    }
    significant_bits_predictor_.Update(context, bit);
    length = (length << 1) + bit;
  }
  *run_length = length;
}

//////////////////////////////////////////////////////////////////
// RunLengthCompressor
//////////////////////////////////////////////////////////////////
namespace{
bool ParseRunLengthCompressorOptions(const std::string& options,
                                     int* parameters) {
  if (options.length() != 13) return false;
  int params[13];
  // CharCoder byte context length
  int byte_context_length = options[0] - '0';
  if (verbosity > 4) {
    std::clog << "Parsing CharCoder options"
              << " byte_context_length=" << byte_context_length
              << std::endl;
  }
  if (byte_context_length < 0 || byte_context_length > 9) return false;
  params[0] = byte_context_length;
  bool ok = true;
  // CharCoder BitHistory option
  ok &= ParseBitHistoryOptions(std::string(options, 1, 3), params + 1);
  // RunLengthCoder BitHistory options
  if (verbosity > 4) {
    std::clog << "Parsing RunlengthCoder options" << std::endl;
  }
  ok &= ParseBitHistoryOptions(std::string(options, 4, 3), params + 4);
  ok &= ParseBitHistoryOptions(std::string(options, 7, 3), params + 7);
  ok &= ParseBitHistoryOptions(std::string(options, 10, 3), params + 10);
  if (ok) memcpy(parameters, params, 13 * sizeof(int));
  return ok;
}
}

bool RunLengthCompressor::Options::Set(const std::string& options) {
  std::string opts = options;
  if (options.length() == 0) opts = "2344255466366";
  if (verbosity > 4) {
    std::clog << "Setting RunLengthCompressor::Options=" << opts << std::endl;
  }
  return ParseRunLengthCompressorOptions(opts, parameters_);
}

std::string RunLengthCompressor::Options::Get() const {
  std::string options;
  for (int i = 0; i < kNumParameters; ++i) {
    assert(parameters_[i] >= 0);
    assert(parameters_[i] <= 9);
    options += static_cast<char>('0' + parameters_[i]);
  }
  return options;
}

void RunLengthCompressor::Write(const char* bytes, size_t n) {
  if (n == 0) return;
  int64 run_length = current_run_length_;
  unsigned char ch = current_char_;
  if (run_length == 0) {
    ch = *bytes++;
    --n;
    run_length = 1;
  }
  for (; n; --n) {
    unsigned char next = *bytes++;
    if (next != ch) {
      EncodeRun(ch, run_length);
      ch = next;
      run_length = 1;
    } else {
      ++run_length;
    }
  }
  current_run_length_ = run_length;
  current_char_ = ch;
}
void RunLengthCompressor::EncodeRun(unsigned char byte, int64 run_length) {
  if (verbosity > 6) {
    std::clog << "Encoding: "
              <<"char=" << int(byte) << " run_length=" << run_length
              << std::endl;
  }
  if (statistics > 0) statistics_collector.RecordRun(byte, run_length);

  char_coder_.Encode(byte);
  run_length_coder_.Encode(byte, run_length);
}

//////////////////////////////////////////////////////////////////
// RunLengthDecompressor
//////////////////////////////////////////////////////////////////
StreamDecompressor* RunLengthDecompressor::Create(const std::string& options) {
  if (options.length() != kNumParameters) return NULL;
  int parameters[kNumParameters];
  if (verbosity > 2) {
    std::clog << "RunLengthDecompressor options=" << options << std::endl;
  }
  bool ok = ParseRunLengthCompressorOptions(options, parameters);
  if (!ok) return NULL;
  return new RunLengthDecompressor(parameters);
}
void RunLengthDecompressor::Read(char* bytes, size_t n) {
  if (n == 0) return;
  int64 run_length = current_run_length_;
  unsigned char ch = current_char_;
  for (; n; --n) {
    if (run_length == 0) {
      DecodeRun(&ch, &run_length);
      assert(run_length > 0);
    }
    *bytes++ = ch;
    --run_length;
  }
  current_run_length_ = run_length;
  current_char_ = ch;
}
void RunLengthDecompressor::DecodeRun(unsigned char* byte, int64* run_length) {
  *byte = char_coder_.Decode();
  *run_length = run_length_coder_.Decode(*byte);
  if (verbosity > 6) {
    std::clog << "Decoded: "
              <<"char=" << int(*byte) << " run_length=" << *run_length
              << std::endl;
  }
}

//////////////////////////////////////////////////////////////////
// RunLengthCompressor::StatisticsCollector
//////////////////////////////////////////////////////////////////
void RunLengthCompressor::StatisticsCollector::StartEncoding() {
  if (statistics <= 0) return;
  memset(num_chars_, 0, 256 * sizeof(int64));
  memset(num_runs_by_char_, 0, 256 * sizeof(int64));
  memset(num_short_runs_by_length_, 0, kLongRunLength * sizeof(int64));
  memset(num_long_runs_by_log_length_, 0,
         (64 - kLogLongRunLength) * sizeof(int64));
}
void RunLengthCompressor::StatisticsCollector::EndEncoding() {
  if (statistics <= 0) return;
  if (statistics > 1) std::clog << "FULL BWT STATISTICS" << std::endl;

  int num_distinct_chars = 0;
  int64 num_runs = 0;
  int64 total_length = 0;
  for (int i = 0; i < 256; ++i) {
      total_length += num_chars_[i];
      num_runs += num_runs_by_char_[i];
      if (num_runs > 0) ++num_distinct_chars;
  }

  for (int i = 0; i < 256; ++i) {
    if (num_chars_[i] > 0) {
      double average_length = (1.0 * num_chars_[i]) / num_runs_by_char_[i];
      double percentage = (100.0 * num_chars_[i]) / total_length;
      if (statistics > 1)
        std::clog << "Char " << i << ": "
                  << num_runs_by_char_[i]
                  << " runs of total length " << num_chars_[i]
                  << std::setprecision(3)
                  << " (" << percentage << "%)"
                  << " and average length " << average_length
                  << std::endl;
    }
  }

  int64 num_short_runs = 0;
  int64 total_length_of_short_runs = 0;
  int64 num_short_runs_by_log_length[kLogLongRunLength];
  memset(num_short_runs_by_log_length, 0, kLogLongRunLength * sizeof(int64));
  for (int length = 1; length < kLongRunLength; ++length) {
    int64 count = num_short_runs_by_length_[length];
    if (count > 0) {
      double percentage = (length * count * 100.0) / total_length;
      if (statistics > 1)
        std::clog << count << " runs of length " << length
                  << std::setprecision(3)
                  << " (" << percentage << "% of total length)"
                  << std::endl;
      num_short_runs += count;
      total_length_of_short_runs += count * length;
      int log_length = 0;
      int len = length;
      while (len >>= 1) ++log_length;
      num_short_runs_by_log_length[log_length] += count;
    }
  }

  for (int log_length = 0; log_length < kLogLongRunLength; ++log_length) {
    int64 count = num_short_runs_by_log_length[log_length];
    if (count > 0) {
      if (statistics > 1)
        std::clog << count << " runs with length in ["
                  << (1 << log_length) << "," << (1 << (log_length + 1)) - 1
                  << "]" << std::endl;
    }
  }

  for (int log_length = kLogLongRunLength; log_length < 64; ++log_length) {
    int64 count = num_long_runs_by_log_length_[log_length - kLogLongRunLength];
    if (count > 0) {
      if (statistics > 1)
        std::clog << count << " runs with length in ["
                  << (1 << log_length) << "," << (1 << (log_length + 1)) - 1
                  << "]" << std::endl;
      int log_length = kLogLongRunLength;
      count >>= kLogLongRunLength;
      while (count >>= 1) ++log_length;
      ++num_short_runs_by_log_length[log_length];
    }
  }

  std::clog << "SUMMARY BWT STATISTICS" << std::endl;
  std::clog << num_distinct_chars << " distinct characters" << std::endl;
  double average_length = (1.0 * total_length) / num_runs;
  std::clog << num_runs << " runs"
            << " of average length " << average_length << std::endl;
  double average_length_of_short_runs =
      (1.0 * total_length_of_short_runs) / num_short_runs;
  std::clog << num_short_runs << " runs shorter than " << kLongRunLength
            << " of total length " << total_length_of_short_runs
            << " and average length " << average_length_of_short_runs
            << std::endl;
  int64 num_long_runs = num_runs - num_short_runs;
  int64 total_length_of_long_runs = total_length - total_length_of_short_runs;
  double average_length_of_long_runs =
      (1.0 * total_length_of_long_runs) / num_long_runs;
  std::clog << num_long_runs << " longer runs"
            << " of total length " << total_length_of_long_runs
            << " and average length " << average_length_of_long_runs
            << std::endl;
}
void RunLengthCompressor::StatisticsCollector::StartBlock() {
  // TODO: Record the statistics at the start of a block.
}
void RunLengthCompressor::StatisticsCollector::EndBlock() {
  // TODO: Compute block statistics by subtracting the statistics
  // at the block start from the current statistics and
  // print them.
}
void RunLengthCompressor::StatisticsCollector::RecordRun(
    unsigned char byte, int64 run_length) {
  if (statistics <= 0) return;
  num_chars_[byte] += run_length;
  ++num_runs_by_char_[byte];
  if (run_length < kLongRunLength) {
    ++num_short_runs_by_length_[run_length];
  } else {
    int log_run_length = kLogLongRunLength;
    run_length >>= kLogLongRunLength;
    while (run_length >>= 1) ++log_run_length;
    ++num_long_runs_by_log_length_[log_run_length - kLogLongRunLength];
  }
}

#endif

}  // namespace dcsbwt
