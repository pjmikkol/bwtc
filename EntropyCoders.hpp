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

namespace bwtc {

/* Base class for all entropy encoders */
class EntropyEncoder {
 public:
  EntropyEncoder() {}
  EntropyEncoder(const std::string& destination, char prob_model) {
    (void) destination;
    (void) prob_model;
  }
  virtual ~EntropyEncoder() {}
  virtual void writeGlobalHeader(char encoding) = 0;
  virtual void encodeData(std::vector<byte>* data, std::vector<uint64>* stats,
                  uint64 data_size) = 0;
  virtual void writeBlockHeader(std::vector<uint64>* stats) = 0;
  virtual void finishBlock(const std::vector<uint32>& LFpowers) = 0;
};

/* Base class for all entropy decoders */
class EntropyDecoder {
 public:
  EntropyDecoder() {}
  EntropyDecoder(const std::string& source) {
    (void) source;
  }
  virtual ~EntropyDecoder() {}
  virtual void readGlobalHeader() = 0;
  virtual std::vector<byte>* decodeBlock(std::vector<uint32>& LFpowers) = 0;
};

EntropyEncoder* giveEntropyEncoder(const std::string& destination, char prob_model);

EntropyDecoder* giveEntropyDecoder(const std::string& source);

} // namespace bwtc

#endif
