/*********************************************************************
 * BlockManager handles the memory usage of a various Block-objects. *
 * With the help of this we can use the same memory multiple times   *
 * without messing up the client code.                               *
 *********************************************************************/

#ifndef BWTC_BLOCK_MANAGER_H_
#define BWTC_BLOCK_MANAGER_H_

#include "block.h"
#include "globaldefs.h"

namespace bwtc {

/* If one wishes to add concurrency this class needs some internal structures
 * for maintaining memory and blocks. */
class BlockManager {
 public:
  BlockManager(uint64 block_size);
  ~BlockManager();
  byte* GetFreeBuffer();
  uint64* GetFreeStats();
  MainBlock* MakeBlock(byte* buffer, uint64* stats, uint64 filled);

 private:
  uint64 block_size_;
  byte* data_buffer_; /* if multiple Block-objects will live at the same *
                       * time, then we need to transform this to array.  */
  uint64* frequency_buffer_; /* above holds also for this */

  BlockManager(const BlockManager&);
  BlockManager& operator=(const BlockManager&);
};

}

#endif
