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

#include <cassert>

#include <vector>

#include "globaldefs.h"

namespace bwtc {

/*******************************************************************
 * MainBlock has two  arrays: block_ and stats_.                   *
 * Both are allocated and deallocated by the client. BlockManager- *
 * class exists for this task.                                     *
 *                                                                 *
 * block_       contains the actual data in MainBlock.             *
 * stats_       contains frequencies for each byte                 *
 *******************************************************************/
class MainBlock {
 public:
  MainBlock(std::vector<byte>* block, std::vector<uint64>* stats,
            uint64 filled);
  ~MainBlock();
  inline uint64 Size() { return filled_; }
  /* begin and end can be used for reading and writing a block.*/
  inline byte* begin()  { return &(*block_)[0]; }
  /* end() returns pointer one past the valid range of its array.
   * Uses filled_ for deducing the value.*/
  inline byte* end() { return &(*block_)[filled_]; }

  /* Append single byte to the end of filled area. */
  inline void Append(byte b) {
    assert(filled_ < block_->size());
    (*block_)[filled_++] = b;
  }
  
  std::vector<byte> *block_;
  std::vector<uint64> *stats_;
  /* Defines range [0, filled_) in array, which will hold the relevant data. */
  uint64 filled_;

 private:
  MainBlock& operator=(const MainBlock& b);
  MainBlock(const MainBlock&);

};



} // namespace bwtc

#endif
