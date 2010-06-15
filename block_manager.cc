#include <vector>

#include "block_manager.h"
#include "block.h"
#include "globaldefs.h"

/* Implementation for trivial BlockManager. */
namespace bwtc {

BlockManager::BlockManager(uint64 block_size, int context_length) :
    block_size_(block_size), data_buffer_(NULL), frequency_buffer_(NULL)
{
  data_buffer_ = new byte[block_size_];
  frequency_buffer_ = new std::vector<uint64>(1 << 8*context_length);
}

BlockManager::~BlockManager() {
  delete [] data_buffer_;
  delete frequency_buffer_;
}

byte* BlockManager::GetFreeBuffer() {
  return data_buffer_;
}

std::vector<uint64>* BlockManager::GetFreeStats() {
  return frequency_buffer_;
}

/* When adding concurrency to the program this one needs re-implementation*/
MainBlock* BlockManager::MakeBlock(byte* data, std::vector<uint64>* stats,
                                   uint64 filled) {
  return new MainBlock(data, stats, filled);
}

} //namespace bwtc
