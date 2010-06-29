#include <cassert>
#include <cstdlib>
#include <ctime>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "../block.h"
#include "../block_manager.h"
#include "../globaldefs.h"
#include "../preprocessors/preprocessor.h"
#include "testdefs.h"

namespace tests {

/*********** begin: TestDefaultPreProcBlockReads() ***********/
void TestDefaultPreProcBlockRead(int fsize, int block_size) {
  /* First write some data into file */
  std::ofstream f(test_fname.c_str());
  std::vector<char> data(fsize, 't');
  std::copy(data.begin(), data.end(), std::ostream_iterator<char>(f));
  f.flush(); f.close();
  /* Then the actual test */
  bwtc::PreProcessor* prepr = bwtc::GivePreProcessor(
      'n', block_size, test_fname);
  bwtc::BlockManager bm(block_size, 1);
  prepr->AddBlockManager(&bm);
  
  int blocks = 0;
  std::streamsize total = 0;
  while (bwtc::MainBlock* b = prepr->ReadBlock()) {
    blocks++;
    total += b->Size();
    delete b;
  }
  delete prepr;
  assert(total == fsize);
  assert( blocks == (fsize / block_size) + 1 ||
          fsize % blocks == 0 ||
         (fsize == 0 && blocks == 0) );
}

void TestDefaultPreProcBlockReads() {
  TestDefaultPreProcBlockRead(1000, 1000);
  srand(time(NULL));
  for (int i = 0; i < 100; ++i) {
    TestDefaultPreProcBlockRead(rand() % 100000, rand() % 100000);
  }
}

/*********** end: TestDefaultPreProcBlockReads() ***********/


} //namespace tests

int main() {
  tests::TestDefaultPreProcBlockReads();
  std::cout << "PreProcessor passed all tests\n";
  return 0;
}
