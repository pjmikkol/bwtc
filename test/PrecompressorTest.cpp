/**
 * @file preproc_test.cc
 * @author Pekka Mikkola <pmikkol@gmail.com>
 *
 * @section LICENSE
 *
 * This file is part of bwtc.
 *
 * bwtc is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * bwtc is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with bwtc.  If not, see <http://www.gnu.org/licenses/>.
 *
 * @section DESCRIPTION
 *
 * Testing common IO capabilites of preprocessor.
 */

#include <cassert>
#include <cstdlib>
#include <ctime>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "../globaldefs.hpp"
#include "../preprocessors/Precompressor.hpp"
#include "../Streams.hpp"
#include "../PrecompressorBlock.hpp"

#undef NDEBUG

namespace bwtc {
int verbosity = 0;
}

namespace tests {

std::string test_fname;

/*********** begin: TestDefaultPreProcBlockReads() ***********/
void TestDefaultPreProcBlockRead(int fsize, int block_size) {
  /* First write some data into file */
  std::ofstream f(test_fname.c_str());
  std::vector<char> data(fsize, 't');
  std::copy(data.begin(), data.end(), std::ostream_iterator<char>(f));
  f.flush(); f.close();
  
  /* Then the actual test */
  bwtc::Precompressor prepr;
  bwtc::RawInStream in(test_fname);

  int blocks = 0;
  std::streamsize total = 0;
  while (true) {
    bwtc::PrecompressorBlock* b = prepr.readBlock(block_size, &in);
    if(b->originalSize() == 0) break;
    blocks++;
    total += b->originalSize();
    delete b;
  }
  assert(total == fsize);
  assert( blocks == (fsize / block_size) + 1 ||
          fsize % blocks == 0 ||
         (fsize == 0 && blocks == 0) );
}

void TestDefaultPreProcBlockReads() {
  TestDefaultPreProcBlockRead(1000, 1000);
  srand(time(0));
  for (int i = 0; i < 100; ++i) {
    TestDefaultPreProcBlockRead(rand() % 100000, rand() % 100000);
  }
}

/*********** end: TestDefaultPreProcBlockReads() ***********/


} //namespace tests

int main(int argc, char **argv) {
  if(argc < 2) {
    std::cout << "Fail: Give test file as a first argument.\n";
    return 1;
  }
  tests::test_fname = argv[1];
  tests::TestDefaultPreProcBlockReads();
  std::cout << "Precompressor passed all tests\n";
  return 0;
}
