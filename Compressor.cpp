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
#include "Streams.hpp"
#include "Profiling.hpp"

#include <string>

namespace bwtc {

Compressor::
Compressor(const std::string& in, const std::string& out,
           const std::string& preprocessing, size_t memLimit, char entropyCoder)
    : m_in(new InStream(in)), m_out(new OutStream(out)),
      m_coder(giveEntropyEncoder(entropyCoder)), m_precompressor(preprocessing),
      m_options(memLimit, entropyCoder) {}

Compressor::
Compressor(InStream* in, OutStream* out,
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
  PROFILE("Compressor::compress");

  size_t compressedSize = 0;
  if(threads != 1) {
    std::cerr << "Supporting only single thread!" << std::endl;
    return 0;
  }
  compressedSize += writeGlobalHeader();

  // More care should be paid for choosing the correct limits
  size_t pbBlockSize = static_cast<size_t>(m_options.memLimit*0.72);
  size_t bwtBlockSize = std::min(static_cast<size_t>(m_options.memLimit*0.18),
                                 static_cast<size_t>(0x80000000 - 2));

  if(m_precompressor.options().size() == 0) pbBlockSize = bwtBlockSize;
  

  size_t preBlocks = 0;
  size_t bwtBlocks = 0;
  while(true) {
    PrecompressorBlock *pb = m_precompressor.readBlock(pbBlockSize, m_in);
    if(pb->originalSize() == 0) {
      delete pb;
      break;
    }
    pb->sliceIntoBlocks(bwtBlockSize);
    ++preBlocks;
    bwtBlocks += pb->slices();

    compressedSize += pb->writeBlockHeader(m_out);

    for(size_t i = 0; i < pb->slices(); ++i) {
      compressedSize += m_coder->
          transformAndEncode(pb->getSlice(i), m_bwtmanager, m_out);
      //TODO: if optimizing overall memory usage now would be time to
      //delete space allocated for i:th slice. However the worst case
      //stays the same
    }


    delete pb;
  }
  
  return compressedSize;
}

} //namespace bwtc
