/* Probably the ugliest hack ever !!*/
#define private public

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
  for(int i = 0; i < size; ++i) {
    data[i] = 0x55;
  }
  data[size-1] = 0x55;
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

void TestWritingAndReadingBlockHeaders() {
  const uint64 header_num = static_cast<uint64>(0x555555555555);

  srand(time(NULL));
  bwtc::Encoder* encoder = new bwtc::Encoder(test_fname, 'n');
  /* Generate data: */
  std::vector<uint64> stats(256);
  for(unsigned i = 0; i < 256; ++i)
    stats[i] = static_cast<uint64>(rand() + 1);
  uint64 header_length;
  std::streampos len_pos = encoder->WriteBlockHeader(&stats, &header_length);
  encoder->out_->Write48bits(header_num, len_pos);
  delete encoder;

  bwtc::Decoder decoder(test_fname);
  std::vector<uint64> context_lengths;
  uint64 compr_len = decoder.ReadBlockHeader(&context_lengths);
  assert(compr_len == header_num);
  assert(context_lengths.size() == 256);
  for(unsigned i = 0; i < 256; ++i) {
    assert(stats[i] == context_lengths[i]);
  }
}

void TestWritingAndReadingHeadersAndSimpleData() {
  const unsigned amount = 30000;

  bwtc::Encoder* encoder = new bwtc::Encoder(test_fname, 'n');
  encoder->WriteGlobalHeader('n','n');
  /* Generate data: */
  std::vector<uint64> stats(256,static_cast<uint64>(0));
  stats[static_cast<byte>('a')] = amount;
  std::vector<byte> data(amount,'a');

  uint64 header_length;
  std::streampos len_pos = encoder->WriteBlockHeader(&stats, &header_length);
  assert(header_length == 5); /* 2 bytes for sentinel and amount*/
  byte* block = &data[0];
  encoder->EncodeRange(block, block + amount);
  encoder->EndContextBlock();
  encoder->Finish();
  encoder->out_->Write48bits(static_cast<uint64>(0xFF), len_pos);
  delete encoder;

  bwtc::Decoder decoder(test_fname);
  char preproc = decoder.ReadGlobalHeader();
  assert(preproc == 'n');
  std::vector<uint64> context_lengths;
  uint64 compr_len = decoder.ReadBlockHeader(&context_lengths);
  assert(compr_len == static_cast<uint64>(0xFF));
  assert(context_lengths.size() == 1);
  decoder.source_->Start();
  for(unsigned i = 0; i < data.size(); ++i) {
    assert(data[i] == decoder.DecodeByte());
  }
  decoder.EndContextBlock();
  int leftover = 0;
  while(decoder.in_->DataLeft()) {
    decoder.in_->ReadByte();
    ++leftover;
  }
  std::cout << leftover << "\n";
}

} //namespace tests

int main() {
  tests::TestArithmeticCoding('n');
  tests::TestWritingAndReadingPackedIntegerList();
  tests::TestPackingIntegers(1000);
  tests::TestBlockArithmeticCoding(1000, 'n');
  tests::TestWritingAndReadingHeadersAndSimpleData();
  tests::TestWritingAndReadingBlockHeaders();
  std::cout << "Encoder and Decoder passed all tests.\n";
}

