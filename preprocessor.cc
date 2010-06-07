#include <iostream>

#include "block.h"
#include "globaldefs.h"
#include "preprocessor.h"
#include "stream.h"

namespace bwtc {

PreProcessor* GivePreProcessor(char choice, int64 block_size) {
  /* Expand this to conform for different PreProcessing algorithms */
  switch (choice) {
    case 'n':
    default:
      return new PreProcessor(block_size);
  }
}

PreProcessor::PreProcessor(int64 block_size) :
    source_(NULL), block_size_(block_size) { }

PreProcessor::~PreProcessor() {}

void PreProcessor::Connect(InStream* source) {
  source_ = source;
}

Block* PreProcessor::ReadBlock() {
  std::vector<char>* block = new std::vector<char>(block_size_);
  /* TODO:
   * streamsize type has as many bits as long. What if user is on 32-bit
   * platform and wants to use too large blocks (amount of bits needed to
   * represent block size is more than 32)?? */
  std::streamsize read = source_->ReadBlock(
      block->begin(), static_cast<std::streamsize>(block_size_) );
  return new Block(block, block->begin() + read);
}

} //namespace bwtc
