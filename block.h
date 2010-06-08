/********************************************************************
 * MainBlock models one block of a data which can have multiple     *
 * context blocks.                                                  *
 *                                                                  *
 * On a compressor pipeline the preprocessor produces MainBlocks    *
 * which are then transformed into cotext blocks by BWT.            *
 *                                                                  *
 * TODO: Do we even need to use MainBlock in decompressing pipeline *
 *       or are byte arrays enough                                  *
 * On a decompressing pipeline ....                                 *
 ********************************************************************/


#ifndef BWTC_BLOCK_H_
#define BWTC_BLOCK_H_

#include <vector>

#include "globaldefs.h"

namespace bwtc {

/*******************************************************************
 * MainBlock has two  arrays: block_ and frequencies_.             *
 * Both are allocated and deallocated by the client. BlockManager- *
 * class exists for this task.                                     *
 *                                                                 *
 * block_       contains the actual data in MainBlock.             *
 * frequencies_ contains frequencies for each byte                 *
 *******************************************************************/
class MainBlock {
 public:
  MainBlock(byte* block, int64* stats, int64 filled);
  ~MainBlock();
  int64 Size() { return filled_; }
  /* begin and end are meant for reading and writing of a block.*/
  byte* begin() { return &block_[0]; }
  /* end() returns pointer one past the valid range of its array.
   * Uses filled_ for deducing the value.*/
  byte* end() { return &block_[filled_]; }

 private:
  byte* block_;
  int64* frequencies_;
  /* Defines range [0, filled_) in array, which will hold the relevant data. */
  int64 filled_;

  MainBlock& operator=(const MainBlock& b);
  MainBlock(const MainBlock&);
  
};

} // namespace bwtc

#endif
