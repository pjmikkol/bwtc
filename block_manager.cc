/**************************************************************************
 *  Copyright 2010, Pekka Mikkola, pjmikkol (at) cs.helsinki.fi           *
 *                                                                        *
 *  This file is part of bwtc.                                            *
 *                                                                        *
 *  bwtc is free software: you can redistribute it and/or modify          *
 *  it under the terms of the GNU General Public License as published by  *
 *  the Free Software Foundation, either version 3 of the License, or     *
 *  (at your option) any later version.                                   *
 *                                                                        *
 *  bwtc is distributed in the hope that it will be useful,               *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *  GNU General Public License for more details.                          *
 *                                                                        *
 *  You should have received a copy of the GNU General Public License     *
 *  along with bwtc.  If not, see <http://www.gnu.org/licenses/>.         *
 **************************************************************************/

#include <vector>

#include "block_manager.h"
#include "block.h"
#include "globaldefs.h"

/* Implementation for trivial BlockManager. */
namespace bwtc {

BlockManager::BlockManager(uint64 block_size, int context_length) :
    block_size_(block_size), data_buffer_(NULL), frequency_buffer_(NULL)
{
  data_buffer_ = new std::vector<byte>(block_size_);
  /* One extra-context for sentinel character */
  frequency_buffer_ = new std::vector<uint64>((1 << 8*context_length) + 1);
}

BlockManager::~BlockManager() {
  delete data_buffer_;
  delete frequency_buffer_;
}

std::vector<byte>* BlockManager::GetFreeBuffer() {
  return data_buffer_;
}

std::vector<uint64>* BlockManager::GetFreeStats() {
  return frequency_buffer_;
}

/* When adding concurrency to the program this one needs re-implementation*/
MainBlock* BlockManager::MakeBlock(std::vector<byte>* data,
                                   std::vector<uint64>* stats, uint64 filled)
{
  return new MainBlock(data, stats, filled);
}

} //namespace bwtc
