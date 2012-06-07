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

namespace bwtc {

PrecompressorBlock::PrecompressorBlock(size_t maxSize, RawInStream* in) {
  m_data.resize(maxSize);
  uint64 read = in->readBlock(&m_data[0], maxSize);
  
}

void PrecompressorBlock::
sliceIntoBlocks(std::vector<BWTBlock>& blocks, uint32 blockSize) {
  size_t begin = 0;
  while(begin < m_used) {
    uint32 bSize = std::min((size_t)blockSize, m_used - begin);
    blocks.push_back(BWTBlock(&m_data[begin], bSize, false));
    begin += bSize;
  }
}

} //namespace bwtc
