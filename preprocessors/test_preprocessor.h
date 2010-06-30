#ifndef BWTC_TEST_PREPROCESSOR_H_
#define BWTC_TEST_PREPROCESSOR_H_

#include <iostream> /* for std::streamsize*/
#include <string>

#include "../block.h"
#include "../block_manager.h"
#include "../globaldefs.h"
#include "../stream.h"
#include "preprocessor.h"

namespace bwtc {
/* This class is meant for testing the different preprocessor-options.
 * It can be done by using the available public-functions such as
 * CompressPairs */
class TestPreProcessor : public PreProcessor {
 public:
  TestPreProcessor(uint64 block_size);
  virtual ~TestPreProcessor();
  /* Reads and preprocesses data to byte array provided by block_manager_*/
  uint64 CompressPairs();
  /* Initialize target-array for reading */
  void InitializeTarget();
  /* Fills the buffer from instream. returns true if something is read*/
  uint64 FillBuffer();

  MainBlock *curr_block_;
};

} // namespace bwtc

#endif
