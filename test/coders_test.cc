#include <cassert>
#include <cstdlib>
#include <ctime>

#include <iostream>
#include <string>

#include "../coders.h"
#include "../globaldefs.h"
#include "../utils.h"
#include "testdefs.h"

namespace bwtc {
int verbosity = 7;
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
  assert('a' == b);
  assert('b' == dec.DecodeByte());
}


void TestBlockArithmeticCoding(int size, char prob_model) {
  std::vector<byte> data(size);
  srand(time(NULL));
  for(int i = 0; i < size; ++i) {
    data[i] = static_cast<byte>(rand() & 0xFF);
  }
  bwtc::Encoder enc(test_fname, prob_model);
  enc.EncodeRange(&data[0], &data[size]);
  enc.Finish();
  bwtc::Decoder dec(test_fname, prob_model);
  dec.Start();
  for(int i = 0; i < size; ++i)
    assert(dec.DecodeByte() == data[i]);
}

void TestPackingIntegers(int times) {
  /* Constants for bit-twiddling */
  static const uint64 kLongOne = 1;

  srand(time(NULL));
  int bytes_needed;
  for(int i = 0; i < times; ++i) {
    uint64 original = static_cast<uint64>(rand());
    uint64 packed = bwtc::PackInteger(original, &bytes_needed);
    uint64 result = bwtc::UnpackInteger(packed);

    assert(result == original);
    assert(original <= (kLongOne << (bytes_needed*8)));
  }
}

} //namespace tests

int main() {
  tests::TestArithmeticCoding('n');
  tests::TestBlockArithmeticCoding(1000, 'n');
  tests::TestPackingIntegers(1000);
  std::cout << "Encoder and Decoder passed all tests.\n";
}

