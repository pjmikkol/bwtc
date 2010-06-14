#include <cassert>

#include <algorithm>
#include <iostream> /* for std::streamsize*/
#include <string>

#include "block_manager.h"
#include "globaldefs.h"
#include "preprocessor.h"
#include "stream.h"

namespace bwtc {

PreProcessor* GivePreProcessor(char choice, uint64 block_size,
                               const std::string& input) {
  PreProcessor* pp;
  /* Expand this to conform for different PreProcessing algorithms */
  switch (choice) {
    case 'n':
    default:
      pp = new PreProcessor(block_size);
  }
  pp->Connect(input);
  return pp;
}

PreProcessor::PreProcessor(uint64 block_size) :
    source_(NULL), block_size_(block_size), block_manager_(NULL) { }

PreProcessor::~PreProcessor() {
  delete source_;
}

void PreProcessor::BuildStats(byte* data, uint64* stats, uint64 size) {
  std::fill(stats, &stats[256], 0);
  for (long i = 0; i < size; ++i) stats[data[i]]++;
}

void PreProcessor::Connect(std::string source_name) {
  source_ = new InStream(source_name);
}

void PreProcessor::AddBlockManager(BlockManager* manager) {
  block_manager_ = manager;
}

MainBlock* PreProcessor::ReadBlock() {
  assert(source_);
  assert(block_manager_);
  byte* to = block_manager_->GetFreeBuffer();
  uint64* stats = block_manager_->GetFreeStats();
  /* TODO:
   * streamsize type has as many bits as long. Since the preprocessor gets
   * blocksize as an uint64 we may end up in problems if user is on 32-bit
   * platform and wants to use too large blocks (amount of bits needed to
   * represent block size is more than 32)?? */
  /*** Stub implementation ***/
  std::streamsize read = source_->ReadBlock(
      to, static_cast<std::streamsize>(block_size_) );
  if (!read) return NULL;
  BuildStats(to, stats, read);
  return block_manager_->MakeBlock(to, stats, static_cast<uint64>(read));
}

} //namespace bwtc
