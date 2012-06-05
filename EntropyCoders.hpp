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

#include "MainBlock.hpp"
#include "globaldefs.hpp"
#include "BitCoders.hpp"
#include "Streams.hpp"
#include "BWTBlock.hpp"
#include "bwtransforms/BWTManager.hpp"
 
namespace bwtc {

/* Base class for all entropy encoders */
class EntropyEncoder {
 public:
  EntropyEncoder() {}
  EntropyEncoder(char prob_model) {
    (void) prob_model;
  }
  virtual ~EntropyEncoder() {}

  virtual size_t transformAndEncode(BWTBlock& block, BWTManager& bwtm,
                                    RawOutStream* out) = 0;

#ifdef ENTROPY_PROFILER
  uint32 m_bytesForCharacters;
  uint32 m_bytesForRuns;
#endif

};

/* Base class for all entropy decoders */
class EntropyDecoder {
 public:
  EntropyDecoder() {}
  EntropyDecoder(const std::string& source) {
    (void) source;
  }
  virtual ~EntropyDecoder() {}
  virtual std::vector<byte>* decodeBlock(std::vector<uint32>& LFpowers) = 0;
};

EntropyEncoder* giveEntropyEncoder(char prob_model);

EntropyDecoder* giveEntropyDecoder(RawInStream* in, char decoder);

} // namespace bwtc

#endif
