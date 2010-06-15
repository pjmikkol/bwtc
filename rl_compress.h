/* Original implementations for BitEncoder and Decoder can be found here: *
 * http://code.google.com/p/dcs-bwt-compressor/                           */

#ifndef DCSBWT_RL_COMPRESS_H__
#define DCSBWT_RL_COMPRESS_H__

#include <string>
#include <cassert>


#include "stream.h"
#include "globaldefs.h" /* Important definitions */

namespace dcsbwt {

using bwtc::verbosity;

/* Probabilities are encoded as integers in [0,kProbabilityScale]
 * with p representing the probability p/kProbabilityScale. */


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
  BitEncoder();
  ~BitEncoder();

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

  /* Measures length of compressed sequence in bytes */
  inline void ResetCounter() { counter_ = 0; }
  inline uint64 Counter() { return counter_; }

 private:
  uint32 low_;
  uint32 high_;
  uint64 counter_;
  bwtc::OutStream* output_;

  inline void EmitByte(unsigned char byte) {
    output_->WriteByte(byte);
    ++counter_;
  }
  BitEncoder(const BitEncoder&);
  BitEncoder& operator=(const BitEncoder&);
};

/*********************************************************************
 * Decompressor for a bit sequence compressed by BitEncoder.         *
 *********************************************************************/
class BitDecoder {
 public:
  BitDecoder();
  ~BitDecoder();

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
  BitDecoder(const BitDecoder&);
  BitDecoder& operator=(const BitDecoder&);
};

}  // namespace dcsbwt

#endif  // DCSBWT_RL_COMPRESS_H__
