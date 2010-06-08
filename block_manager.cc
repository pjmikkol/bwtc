#include "block_manager.h"
#include "block.h"
#include "globaldefs.h"

/* Implementation for trivial BlockManager. */
namespace bwtc {

BlockManager::BlockManager(int64 block_size) :
    block_size_(block_size), data_buffer_(NULL), frequency_buffer_(NULL)
{
  data_buffer_ = new byte[block_size_];
  frequency_buffer_ = new int64[256];
}

BlockManager::~BlockManager() {
  delete [] data_buffer_;
  delete [] frequency_buffer_;
}

byte* BlockManager::GetFreeBuffer() {
  return data_buffer_;
}

int64* BlockManager::GetFreeStats() {
  return frequency_buffer_;
}

/* When adding concurrency to the program this one needs re-implementation*/
MainBlock* BlockManager::MakeBlock(byte* data, int64* stats, int64 filled) {
  return new MainBlock(data, stats, filled);
}

} //namespace bwtc
