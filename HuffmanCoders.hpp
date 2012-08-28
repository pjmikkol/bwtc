/**
 * @file HuffmanCoders.hpp
 * @author Dominik Kempa <dominik.kempa@cs.helsinki.fi>
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
 * Header for simple Huffman coder.
 *
 */

#ifndef BWTC_HUFFMAN_CODERS_HPP_
#define BWTC_HUFFMAN_CODERS_HPP_

#include "EntropyCoders.hpp"
#include "globaldefs.hpp"
#include "BitCoders.hpp"
#include "Streams.hpp"
#include "BWTBlock.hpp"

#include <iostream>
#include <string>
#include <vector>


namespace bwtc {

class HuffmanEncoder : public EntropyEncoder {
 public:
  HuffmanEncoder();
  ~HuffmanEncoder();

  size_t transformAndEncode(BWTBlock& block, BWTManager& bwtm,
                            OutStream* out);
  
  void encodeData(const byte* data, const std::vector<uint32>& stats,
                  uint32 blockSize, OutStream* out);
  void writeBlockHeader(const BWTBlock& b, std::vector<uint32>& stats,
                        OutStream* out);
  void finishBlock(OutStream* out);

  void writePackedInteger(uint64 packed_integer, OutStream* out);


 private:
  long int m_headerPosition;
  uint64 m_compressedBlockLength;

  void serializeShape(uint32 *clen, std::vector<bool> &vec);
  HuffmanEncoder(const HuffmanEncoder&);
  HuffmanEncoder& operator=(const HuffmanEncoder&);
};

class HuffmanDecoder : public EntropyDecoder {
 public:
  HuffmanDecoder();
  ~HuffmanDecoder();

  uint64 readPackedInteger(InStream* in);
  void decodeBlock(BWTBlock& block, InStream* in);
  uint64 readBlockHeader(BWTBlock& block, std::vector<uint64>* stats,
                         InStream* in);

 private:

  size_t deserializeShape(InStream &input, uint32 *clen);
  HuffmanDecoder(const HuffmanDecoder&);
  HuffmanDecoder& operator=(const HuffmanDecoder&);
};

} // namespace bwtc

#endif
