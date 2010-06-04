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
  
}

} //namespace bwtc
