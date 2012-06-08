/**
 * @file PrecompressorBlock.cpp
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
 * Implementation of Precompressor-block.
 */

#include "globaldefs.hpp"
#include "PrecompressorBlock.hpp"
#include "Streams.hpp"

#include <cassert>

namespace bwtc {

PrecompressorBlock::PrecompressorBlock(size_t maxSize, RawInStream* in) {
  m_data.resize(maxSize+1);
  uint64 read = in->readBlock(&m_data[0], maxSize);
  m_data.resize(read+1);
  m_originalSize = m_used = read;
}

void PrecompressorBlock::setSize(size_t size) {
  assert(size > 0);
  m_used = size;
  m_data.resize(size+1);
}

size_t PrecompressorBlock::writeBlockHeader(RawOutStream* out) const {
  
}

void PrecompressorBlock::sliceIntoBlocks(size_t blockSize) {
  //Have at least one additional byte for the sentinel of BWT
  assert(m_used < m_data.size());
  assert(blockSize < (0x80000000 - 2));
  m_bwtBlocks.clear();
  size_t begin = 0;
  while(begin < m_used) {
    uint32 bSize = std::min(blockSize, m_used - begin);
    m_bwtBlocks.push_back(BWTBlock(&m_data[begin], bSize, false));
    begin += bSize;
  }
}


BWTBlock& PrecompressorBlock::getSlice(int i) {
  assert(i >= 0);
  assert((size_t)i < m_bwtBlocks.size());
  return m_bwtBlocks[i];
}

} //namespace bwtc
