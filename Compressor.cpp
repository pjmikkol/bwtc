/**
 * @file Compressor.cpp
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
 * Implementation of Compressor.
 */

#include "Compressor.hpp"

#include <string>

namespace bwtc {

Compressor::
Compressor(const std::string& in, const std::string& out, size_t memLimit,
           char entropyCoder)
    : m_in(new RawInStream(in)), m_out(new RawOutStream(out)),
      m_coder(giveEntropyEncoder(entropyCoder)),
      m_preprocessor(0), m_options(memLimit, entropyCoder) {}

Compressor::Compressor(RawInStream* in, RawOutStream* out, size_t memLimit,
                       char entropyCoder)
    : m_in(in), m_out(out), m_coder(giveEntropyEncoder(entropyCoder)),
      m_preprocessor(0), m_options(memLimit, entropyCoder) {}

Compressor::~Compressor() {
  delete m_in;
  delete m_out;
  delete m_coder;
  delete m_preprocessor;
}

size_t Compressor::writeGlobalHeader() {
  m_out->writeByte(static_cast<byte>(m_options.entropyCoder));
  return 1;
}

void Compressor::setPreprocessor(const std::string& parameters) {
  delete m_preprocessor;
  m_preprocessor = new Preprocessor(parameters);
}

size_t Compressor::compress(size_t threads) {
  size_t compressedSize = 0;
  if(threads != 1) {
    std::cerr << "Supporting only single thread!" << std::endl;
    return 0;
  }
  compressedSize += writeGlobalHeader();

  // More care should be paid for choosing the correct limits
  size_t pbBlockSize = static_cast<size_t>(m_options.memLimit*0.72);
  size_t bwtBlockSize = static_cast<size_t>(memLimit*0.18);
  //std::vector<byte> temp(m_options.memLimit*0.24); for preprocessor

  while(true) {
    PrecompressorBlock *pb = m_preprocessor.readBlock(pbBlockSize);
    if(pb->originalSize() == 0) {
      delete pb;
      break;
    }
    //Write pb-header: PrecompressorBlock or Compressor ?

    std::vector<BWTBlock> bwtBlocks;
    size_t begin = 0;
    while(begin < pb->size()) {
      bwtBlocks.push_back(pb->makeBlock(begin, bwtBlockSize));
      begin += bwtBlocks.back().size();
    }

    for(size_t i = 0; i < bwtBlocks.size(); ++i) {
      compressedSize += m_coder->
          transformAndEncode(bwtBlocks[i], m_bwtManager, m_out);
    }


    delete pb;
  }
  
  
  return compressedSize;
}

} //namespace bwtc
