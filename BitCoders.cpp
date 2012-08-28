/**
 * @file BitCoders.cpp
 * @author http://code.google.com/p/dcs-bwt-compressor/ Original implementations 
 * @author Pekka Mikkola <pmikkol@gmail.com>
 *
 * @section LICENSE
 *
 * This file is part of bwtc.
 *
 * bwtc is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * bwtc is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with bwtc.  If not, see <http://www.gnu.org/licenses/>.
 *
 * @section DESCRIPTION
 *
 * Implementations for BitEncoder and BitDecoder.
 */

#include <cassert>

#include <string>

#include "BitCoders.hpp"
#include "Streams.hpp"

namespace dcsbwt {

extern int statistics;

/**
 * BitEncoder and BitDecoder maintain a range [m_low,m_high] (of uint32's).
 * The next four bytes of the compressed output/input represent
 * some value in this range.
 * For each bit (en/de)coded, the range is split into two parts
 * with the sizes proportional to the probabilities of 1-bit and 0-bit.
 * BitEncoder chooses one the two subranges as the new range [m_low,m_high].
 * based on the actual value of the bit.
 * BitDecoder determines the bit by finding out, which subrange
 * contains the value represented by the next four input bytes.
 * Whenever m_low and m_high share the most significant byte,
 * the range is expanded by the factor 256 shifting the shared byte out.
 * BitEncoder emits the shared byte. BitDecoder shifts the shared byte
 * out of the next four input bytes that it is holding and reads in
 * a new input byte.
 */

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
    : m_low(0), m_high(0xFFFFFFFF), m_counter(0), m_output(NULL) {}

BitEncoder::~BitEncoder() { }

void BitEncoder::encode(bool bit, Probability probability_of_one) {
  if (verbosity > 7) {
    std::clog << "Encoding bit " << int(bit) << " with probability "
              << (bit ? probability_of_one :
                  kProbabilityScale - probability_of_one)
              << std::endl;
  }
  uint32 split = Split(m_low, m_high, probability_of_one);
  /* Choose the subrange. */
  if (bit) m_high = split; else m_low = split + 1;
  while (((m_low ^ m_high) & 0xFF000000) == 0) {
    emitByte(m_low >> 24);
    m_low <<= 8;
    m_high = (m_high << 8) + 255;
  }
  assert(m_low < m_high);
}

void BitEncoder::finish() {
  /* Emit 4 bytes representing any value in [m_low,m_high]
   * These are needed for BitDecoder lookahead. */
  emitByte(m_low >> 24);
  emitByte(255);
  emitByte(255);
  emitByte(255);
  m_output->flush();
  /* Prepare to encode another sequence. */
  m_low = 0;
  m_high = 0xFFFFFFFF;
}

BitDecoder::BitDecoder() :
    m_low(0), m_high(0xFFFFFFFF), m_next(0), m_input(NULL) {}

BitDecoder::~BitDecoder() { }

void BitDecoder::start() {
  m_low = 0;
  m_high = 0xFFFFFFFF;
  m_next = readByte();
  m_next = (m_next << 8) + readByte();
  m_next = (m_next << 8) + readByte();
  m_next = (m_next << 8) + readByte();
}

bool BitDecoder::decode(Probability probability_of_one) {
  uint32 split = Split(m_low, m_high, probability_of_one);
  bool bit = (m_next <= split);
  if (bit) m_high = split; else m_low = split + 1;
  while (((m_low ^ m_high) & 0xFF000000) == 0) {
    m_low <<= 8;
    m_high = (m_high << 8) + 255;
    m_next = (m_next << 8) + readByte();
  }
  assert(m_low < m_high);
  assert(m_next >= m_low);
  assert(m_next <= m_high);
  if (verbosity > 7) {
    std::clog << "Decoded bit " << int(bit) << " with probability "
              << (bit ? probability_of_one :
                  kProbabilityScale - probability_of_one)
              << std::endl;
  }
  return bit;
}

}  // namespace dcsbwt
