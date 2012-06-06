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
#include "PrecompressorBlock.hpp"

#include <string>

namespace bwtc {

Compressor::
Compressor(const std::string& in, const std::string& out,
           const std::string& preprocessing, size_t memLimit, char entropyCoder)
    : m_in(new RawInStream(in)), m_out(new RawOutStream(out)),
      m_coder(giveEntropyEncoder(entropyCoder)), m_precompressor(preprocessing),
      m_options(memLimit, entropyCoder) {}

Compressor::
Compressor(RawInStream* in, RawOutStream* out,
           const std::string& preprocessing, size_t memLimit, char entropyCoder)
    : m_in(in), m_out(out), m_coder(giveEntropyEncoder(entropyCoder)),
      m_precompressor(preprocessing), m_options(memLimit, entropyCoder) {}

Compressor::~Compressor() {
  delete m_in;
  delete m_out;
  delete m_coder;
}

size_t Compressor::writeGlobalHeader() {
  m_out->writeByte(static_cast<byte>(m_options.entropyCoder));
  return 1;
}

void Compressor::initializeBwtAlgorithm(char choice, uint32 startingPoints) {
  m_bwtmanager.initialize(choice);
  m_bwtmanager.setStartingPoints(startingPoints);
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
  size_t bwtBlockSize = static_cast<size_t>(m_options.memLimit*0.18);
  //std::vector<byte> temp(m_options.memLimit*0.24); for preprocessor

  while(true) {
    PrecompressorBlock *pb = m_precompressor.readBlock(pbBlockSize, m_in);
    if(pb->originalSize() == 0) {
      delete pb;
      break;
    }
    //Write pb-header: PrecompressorBlock or Compressor ?

    std::vector<BWTBlock> bwtBlocks;
    pb->sliceIntoBlocks(bwtBlocks, bwtBlockSize);

    for(size_t i = 0; i < bwtBlocks.size(); ++i) {
      compressedSize += m_coder->
          transformAndEncode(bwtBlocks[i], m_bwtmanager, m_out);
      //TODO: if optimizing memory usage now would be time to
      //delete space allocated for bwtBlocks[i]
    }


    delete pb;
  }
  
  
  return compressedSize;
}

} //namespace bwtc
