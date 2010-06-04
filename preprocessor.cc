#include "preprocessor.h"
#include "stream.h"

namespace bwtc {

PreProcessor* GivePreProcessor(char choice) {
  /* Expand this to conform for different PreProcessing algorithms */
  switch (choice) {
    case 'n':
    default:
      return new PreProcessor();
  }
}


PreProcessor::PreProcessor() :
    source_(NULL) { }

PreProcessor::~PreProcessor() {}

void PreProcessor::Connect(InStream* source) {
  source_ = source;
}


} //namespace bwtc
