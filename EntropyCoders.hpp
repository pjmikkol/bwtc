/**
 * @file EntropCoders.hpp
 * @author Dominik Kempa <dominik.kempa@cs.helsinki.fi>
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
 * Header for base classes for entropy coders.
 *
 */

#ifndef BWTC_ENTROPY_CODERS_HPP_
#define BWTC_ENTROPY_CODERS_HPP_

#include <iostream>
#include <string>
#include <vector>

#include "globaldefs.hpp"
#include "BitCoders.hpp"
#include "Streams.hpp"
#include "BWTBlock.hpp"
#include "bwtransforms/BWTManager.hpp"
 
namespace bwtc {

/* Base class for all entropy encoders */
class EntropyEncoder {
 public:
  virtual ~EntropyEncoder() {}

  /**Entropy coder should call BWTManager to calculate transform and after that
   * compress the transformed string. By calling BWTManager, encoder can
   * have some additional information about the string which can be calculated
   * during the computation of BWT. With this one pass over the transformed string
   * can be avoided.
   *
   * @param block Block to transform and encode.
   * @param bwtm BWTManager to make a transform.
   * @param out Stream where the encoded result is written.
   * @return total bytes for the compressed block (including header)
   */
  virtual size_t transformAndEncode(BWTBlock& block, BWTManager& bwtm,
                                    OutStream* out) = 0;

#ifdef ENTROPY_PROFILER
  uint32 m_bytesForCharacters;
  uint32 m_bytesForRuns;
#endif

};

/* Base class for all entropy decoders */
class EntropyDecoder {
 public:
  virtual ~EntropyDecoder() {}
  virtual void decodeBlock(BWTBlock& block, InStream* in) = 0;
};

EntropyEncoder* giveEntropyEncoder(char encoder);

EntropyDecoder* giveEntropyDecoder(char decoder);

} // namespace bwtc

#endif
