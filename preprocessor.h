#ifndef BWTC_PREPROCESSOR_H_
#define BWTC_PREPROCESSOR_H_

#include "stream.h"

namespace bwtc {

//TODO: Make this class to abstract class (interface for real implementations)
class PreProcessor {
 public:
  PreProcessor();
  ~PreProcessor();
  void Connect(InStream* source);

 private:
  InStream* source_;

  PreProcessor& operator=(const PreProcessor& p);
  PreProcessor(const PreProcessor&);
};

/* This function returns chosen preprocessor */ 
PreProcessor* GivePreProcessor(char choice);

} // namespace bwtc

#endif
