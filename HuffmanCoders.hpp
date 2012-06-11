/**
 * @file HuffmanCoders.hpp
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
 * Header for simple Huffman coder.
 *
 */

#ifndef BWTC_HUFFMAN_CODERS_HPP_
#define BWTC_HUFFMAN_CODERS_HPP_

#include <iostream>
#include <string>
#include <vector>

#include "EntropyCoders.hpp"
#include "globaldefs.hpp"
#include "BitCoders.hpp"
#include "Streams.hpp"

namespace bwtc {

class HuffmanEncoder : public EntropyEncoder {
 public:
  HuffmanEncoder(const std::string& destination, char prob_model);
  ~HuffmanEncoder();
  void writeGlobalHeader(char encoding);
  void encodeData(const byte* data, std::vector<uint64>* stats,
                  uint64 data_size);
  void writeBlockHeader(std::vector<uint64>* stats);

  void writePackedInteger(uint64 packed_integer);
  int writeTrailer(const std::vector<uint32>& LFpowers);
  void finishBlock(const std::vector<uint32>& LFpowers);

 private:
  RawOutStream* m_out;
  std::streampos m_headerPosition;
  uint64 m_compressedBlockLength;

  void serializeShape(uint32 *clen, std::vector<bool> &vec);
  HuffmanEncoder(const HuffmanEncoder&);
  HuffmanEncoder& operator=(const HuffmanEncoder&);
};

class HuffmanDecoder : public EntropyDecoder {
 public:
  HuffmanDecoder(const std::string& source);
  HuffmanDecoder(RawInStream *in);
  ~HuffmanDecoder();
  void readGlobalHeader();
  uint64 readPackedInteger();
  std::vector<byte>* decodeBlock(std::vector<uint32>& LFpowers);
  uint64 readBlockHeader(std::vector<uint64>* stats);

 private:
  RawInStream* m_in;

  size_t deserializeShape(RawInStream &input, uint32 *clen);
  HuffmanDecoder(const HuffmanDecoder&);
  HuffmanDecoder& operator=(const HuffmanDecoder&);
};

} // namespace bwtc

#endif
