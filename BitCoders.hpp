/**
 * @file BitCoders.hpp
 * @author http://code.google.com/p/dcs-bwt-compressor/ Original implementations 
 * @author Pekka Mikkola <pjmikkol@cs.helsinki.fi>
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
 * Header for BitEncoder and BitDecoder.
 */

#ifndef DCSBWT_RL_COMPRESS_HPP__
#define DCSBWT_RL_COMPRESS_HPP__

#include <string>
#include <cassert>


#include "Streams.hpp"
#include "globaldefs.hpp" /* Important definitions */

using bwtc::uint64;
using bwtc::uint32;
using bwtc::byte;
using bwtc::Probability;
using bwtc::kProbabilityScale;
using bwtc::kLogProbabilityScale;

namespace dcsbwt {

using bwtc::verbosity;

/* Probabilities are encoded as integers in [0,kProbabilityScale]
 * with p representing the probability p/kProbabilityScale. */

/**
 * Entropy compressor for a sequence of bits.
 * Given a probility distribution for each bit in the sequence,
 * ouputs a compressed representation of the sequence.
 * Uses range encoding (arithmetic encoding).
 * The expected length of the compressed representation is
 * close to the information theoretic lower bound.
 */

class BitEncoder {
 public:
  BitEncoder();
  ~BitEncoder();

  void connect(bwtc::RawOutStream* out) { m_output = out; }
  //TODO: Figure out what Disconnect should do if needed
  //bwtc::OutStream* Disconnect() { return output_.Disconnect(); }

  /* Append a bit to the sequence.
   * Even input cases bit==1(true), probability_of_one==0 and
   * bit==0(false), probability_of_one==kProbabilityScale are legal.
   * In these cases, up to four bytes of output is generated
   * from the single bit.  */
  void encode(bool bit, Probability probability_of_one);

  /* Must be called to finish the encoding of a sequence.
   * After the call, BitEncoder is ready to start encoding a new sequence. */
  void finish();

  /* Measures length of compressed sequence in bytes */
  inline void resetCounter() { m_counter = 0; }
  inline uint64 counter() { return m_counter; }

 private:
  uint32 m_low;
  uint32 m_high;
  uint64 m_counter;
  bwtc::RawOutStream* m_output;

  inline void emitByte(unsigned char byte) {
    m_output->writeByte(byte);
    ++m_counter;
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
  void connect(bwtc::RawInStream* in) { m_input = in; }
  //TODO: Do we need Disconnect()?
  //bwtc::InStream* Disconnect() { return input_.Disconnect(); }

  /* start() must be called to start the decoding of a sequence.
   * Nothing needs to be called to finish the decoding, but start()
   * must be called again when starting to decode a new sequence.  */
  void start();

  /* Get the next bit of the uncompressed sequence.
   * The probability distribution must be the same as the one used
   * when encoding the same bit. */
  bool decode(Probability probability_of_one);

 private:
  uint32 m_low;
  uint32 m_high;
  uint32 m_next;
  bwtc::RawInStream* m_input;

  byte readByte() { return m_input->readByte(); }
  BitDecoder(const BitDecoder&);
  BitDecoder& operator=(const BitDecoder&);
};

}  // namespace dcsbwt

#endif  // DCSBWT_RL_COMPRESS_H__
