#include <cassert>
#include <cstdlib>
#include <ctime>

#include <iostream>
#include <string>
#include <vector>

#include "../coders.h"
#include "../globaldefs.h"
#include "testdefs.h"

namespace bwtc {
int verbosity = 9;
}

namespace tests {


void TestArithmeticCoding(char prob_model) {
  bwtc::Encoder enc(test_fname, prob_model);
  enc.EncodeByte('a');
  enc.EncodeByte('b');
  enc.Finish();
  bwtc::Decoder dec(test_fname, prob_model);
  dec.Start();
  byte b = dec.DecodeByte();
  std::cout << int(b) << "\n";
  assert(b == dec.DecodeByte());
  assert('b' == dec.DecodeByte());
}


void TestSimpleArithmeticCoding(int size, char prob_model) {
  std::vector<byte> data(size);
  srand(time(NULL));
  for(int i = 0; i < size; ++i) {
    data[i] = static_cast<byte>(rand() % 0xFF);
  }
  bwtc::Encoder enc(test_fname, prob_model);
  enc.EncodeBlock(&data[0], &data[size]);
  enc.Finish();
  bwtc::Decoder dec(test_fname, prob_model);
  dec.Start();
  for(int i = 0; i < size; ++i) {
    assert(dec.DecodeByte() == data[i]);
  }
}

} //namespace tests

int main() {
  tests::TestArithmeticCoding('n');
  std::cout << "Encoder and Decoder passed all tests.\n";
}

