#include <cassert>
#include <cstdlib>
#include <ctime>

#include <iostream>
#include <string>
#include <vector>

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

void TestWritingAndReadingPackedIntegerList() {
  static const uint64 kErrorMask = static_cast<uint64>(1) << 63;

  srand(time(NULL));
  int symbols = rand() & 0x00FFFF;
  bwtc::Encoder* enc = new bwtc::Encoder(test_fname,'n');
  std::vector<uint64> data(symbols);
  int bytes;
  for(unsigned i = 0; i < data.size(); ++i) {
    data[i] = static_cast<uint64>(rand());
    uint64 packed_integer = bwtc::PackInteger(data[i], &bytes);
    enc->WritePackedInteger(packed_integer, bytes);
  }
  enc->FinishBlockHeader();
  delete enc;
  bwtc::Decoder dec(test_fname,'n');
  unsigned i = 0;
  while (1) {
    uint64 value = dec.ReadPackedInteger();
    if(value & kErrorMask) break;
    assert(i < data.size());
    assert(data[i] == bwtc::UnpackInteger(value));
    ++i;
  }
  assert(i == data.size());
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
  tests::TestWritingAndReadingPackedIntegerList();
  std::cout << "Encoder and Decoder passed all tests.\n";
}

