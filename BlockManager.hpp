/**
 * @file BlockManager.hpp
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
 * Header file for BlockManager.
 */


#ifndef BWTC_BLOCK_MANAGER_H_
#define BWTC_BLOCK_MANAGER_H_

#include "MainBlock.hpp"
#include "globaldefs.hpp"

namespace bwtc {

/* If one wishes to add concurrency this class needs some internal structures
 * for maintaining memory and blocks. TODO: Interface needs some rewriting */


/**
 * BlockManager handles the memory usage of a various MainBlock-objects.
 * With the help of this we can use the same memory multiple times
 * without worries of the memory management in code using blocks.
 */
class BlockManager {
 public:
  BlockManager(uint64 block_size, int context_length);
  ~BlockManager();
  std::vector<byte>* getFreeBuffer();
  std::vector<uint64>* getFreeStats();
  MainBlock* makeBlock(std::vector<byte>* buffer, std::vector<uint64>* stats,
                       uint64 filled);

 private:
  uint64 m_blockSize;
  /* If multiple Block-objects will exist  at the same time, then we need to
   * transform this to array.  */
  std::vector<byte>* m_dataBuffer; 
  /* above holds also for this */
  std::vector<uint64>* m_frequencyBuffer; 

  BlockManager(const BlockManager&);
  BlockManager& operator=(const BlockManager&);
};

}

#endif
