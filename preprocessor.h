#ifndef BWTC_PREPROCESSOR_H_
#define BWTC_PREPROCESSOR_H_

#include <string>

#include "block.h"
#include "globaldefs.h"
#include "stream.h"

namespace bwtc {

//TODO: Make this class to abstract class (interface for real implementations)
class PreProcessor {
 public:
  PreProcessor(int64 block_size);
  ~PreProcessor();
  void Connect(std::string source_name);
  Block* ReadBlock();

 private:
  InStream* source_;
  int64 block_size_;

  PreProcessor& operator=(const PreProcessor& p);
  PreProcessor(const PreProcessor&);
};

/* This function returns chosen preprocessor */ 
PreProcessor* GivePreProcessor(char choice, int64 block_size);

} // namespace bwtc

#endif
