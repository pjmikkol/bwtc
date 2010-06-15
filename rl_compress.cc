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

}  // namespace dcsbwt
