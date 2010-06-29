#include <cassert>

#include <vector>

#include "../block.h"
#include "../block_manager.h"
#include "../globaldefs.h"
#include "../stream.h"
#include "preprocessor.h"
#include "test_preprocessor.h"

namespace bwtc {

TestPreProcessor::TestPreProcessor(uint64 block_size) :
    PreProcessor(block_size), curr_block_(NULL) {}

TestPreProcessor::~TestPreProcessor() {}

uint64 TestPreProcessor::CompressPairs() {
  uint64 filled = CompressCommonPairs(&(*curr_block_->block_)[0],
                                      curr_block_->filled_);
  uint64 result = curr_block_->filled_ - filled;
  curr_block_->filled_ = filled;
  return result;
}

void TestPreProcessor::InitializeTarget() {
  assert(block_manager_);
  std::vector<byte>* target = block_manager_->GetFreeBuffer();
  std::vector<uint64>* stats = block_manager_->GetFreeStats();
  curr_block_ = block_manager_->MakeBlock(target, stats, 0UL);
}

uint64 TestPreProcessor::FillBuffer() {
  // TODO: Check the types
  assert(source_);
  assert(curr_block_);
  assert(block_size_ == curr_block_->block_->size());
  if (curr_block_->filled_ == block_size_) return 0;
  std::streamsize read = source_->ReadBlock(
      curr_block_->begin() + curr_block_->Size(),
      static_cast<std::streamsize>(block_size_ - curr_block_->Size()));
  curr_block_->filled_ += read;
  return read;
}


} //namespace bwtc
