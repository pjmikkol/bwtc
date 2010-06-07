#include <cassert>
#include <cstdlib>
#include <ctime>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "../block.h"
#include "../globaldefs.h"
#include "../preprocessor.h"
#include "testdefs.h"

namespace tests {

/*********** begin: TestDefaultPreProcBlockReads() ***********/
void TestDefaultPreProcBlockRead(int fsize, int block_size) {
  std::ofstream f(test_fname.c_str());
  std::vector<char> data(fsize, 't');
  std::copy(data.begin(), data.end(), std::ostream_iterator<char>(f));
  f.flush(); f.close();
  bwtc::PreProcessor* prepr = bwtc::GivePreProcessor('n', block_size);
  prepr->Connect(test_fname);
  int blocks = 0;
  while (bwtc::Block* b = prepr->ReadBlock()) {
    blocks++;
    delete b;
  }
  delete prepr;
  assert(blocks*block_size >= fsize);
  assert(blocks == (fsize / block_size) + 1 || !fsize);
}

void TestDefaultPreProcBlockReads() {
  TestDefaultPreProcBlockRead(0, 1000);
  srand(time(NULL));
  for (int i = 0; i < 10; ++i) {
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
