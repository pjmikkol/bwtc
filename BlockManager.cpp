/**
 * @file BlockManager.cpp
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
 * Implementation for trivial BlockManager.
 */

#include <vector>

#include "BlockManager.hpp"
#include "MainBlock.hpp"
#include "globaldefs.hpp"

namespace bwtc {

BlockManager::BlockManager(uint64 block_size, int context_length) :
    m_blockSize(block_size),
    m_dataBuffer(new std::vector<byte>(m_blockSize)),
    // One extra-context for sentinel character
    m_frequencyBuffer(new std::vector<uint64>((1 << 8*context_length) + 1))
{
}

BlockManager::~BlockManager() {
  delete m_dataBuffer;
  delete m_frequencyBuffer;
}

std::vector<byte>* BlockManager::GetFreeBuffer() {
  return m_dataBuffer;
}

std::vector<uint64>* BlockManager::GetFreeStats() {
  return m_frequencyBuffer;
}

/* When adding concurrency to the program this one needs re-implementation*/
MainBlock* BlockManager::MakeBlock(std::vector<byte>* data,
                                   std::vector<uint64>* stats, uint64 filled)
{
  return new MainBlock(data, stats, filled);
}

} //namespace bwtc
