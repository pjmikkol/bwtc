/**
 * @file coders_test.cc
 * @author Pekka Mikkola <pjmikkol@cs.helsinki.fi>
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
 * Testing of arithemic coder and decoder.
 *
 */

/* Probably the ugliest hack ever !!*/
#define private public

#include <cassert>
#include <cstdlib>
#include <ctime>

#include <iostream>
#include <string>
#include <vector>

#include "../Coders.hpp"
#include "../globaldefs.hpp"
#include "../Utils.hpp"

namespace bwtc {
int verbosity = 1;
}

namespace tests {

std::string test_fname;

void TestArithmeticCoding(char prob_model) {
  bwtc::Encoder enc(test_fname, prob_model);
  enc.encodeByte('a');
  enc.encodeByte('b');
  enc.m_destination->finish();
  bwtc::Decoder dec(test_fname, prob_model);
  dec.start();
  byte b = dec.decodeByte();
  assert('a' == b);
  assert('b' == dec.decodeByte());
}

void TestWritingAndReadingPackedIntegerList() {
  static const uint64 kErrorMask = static_cast<uint64>(1) << 63;

  srand(time(NULL));
  int symbols = rand() & 0x00FFFF;
  bwtc::Encoder* enc = new bwtc::Encoder(test_fname,'n');
  std::vector<uint64> data(symbols);
  data.push_back(0xF0000FF0);
  int bytes;
  for(unsigned i = 0; i < data.size(); ++i) {
    data[i] = static_cast<uint64>(rand());
    uint64 packed_integer = utils::packInteger(data[i], &bytes);
    enc->writePackedInteger(packed_integer);
  }
  enc->finishBlockHeader();
  delete enc;
  bwtc::Decoder dec(test_fname,'n');
  unsigned i = 0;
  while (1) {
    uint64 value = dec.readPackedInteger();
    if(value & kErrorMask) break;
    assert(i < data.size());
    assert(data[i] == utils::unpackInteger(value));
    ++i;
  }
  assert(i == data.size());
}

#if 0
void TestBlockArithmeticCoding(int size, char prob_model) {
  std::vector<byte> data(size);
  for(int i = 0; i < size; ++i) {
    data[i] = 0x55;
  }
  data[size-1] = 0x55;
  bwtc::Encoder enc(test_fname, prob_model);
  enc.encodeRange(&data[0], &data[size]);
  enc.FinishBlock(0x0F);
  bwtc::Decoder dec(test_fname, prob_model);
  dec.Start();
  for(int i = 0; i < size; ++i)
    assert(dec.decodeByte() == data[i]);
}

void TestPackingIntegers(int times) {
  /* Constants for bit-twiddling */
  static const uint64 kLongOne = 1;

  srand(time(NULL));
  int bytes_needed;
  for(int i = 0; i < times; ++i) {
    uint64 original = static_cast<uint64>(rand());
    uint64 packed = utils::packInteger(original, &bytes_needed);
    uint64 result = utils::unpackInteger(packed);

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
  encoder->writeBlockHeader(&stats);
  encoder->out_->write48bits(header_num, 0x0F0F);
  delete encoder;

  bwtc::Decoder decoder(test_fname);
  std::vector<uint64> context_lengths;
  uint64 compr_len = decoder.readBlockHeader(&context_lengths);
  assert(compr_len == header_num);
  assert(context_lengths.size() == 256);
  for(unsigned i = 0; i < 256; ++i) {
    assert(stats[i] == context_lengths[i]);
  }
}


void TestWritingAndReadingHeadersAndSimpleData() {
  const unsigned amount = 300;

  bwtc::Encoder* encoder = new bwtc::Encoder(test_fname, 'n');
  encoder->writeGlobalHeader('n','n');
  /* Generate data: */
  std::vector<uint64> stats(256,static_cast<uint64>(0));
  stats[static_cast<byte>('a')] = amount;
  std::vector<byte> data(amount,'a');

  uint64 header_length;
  std::streampos len_pos = encoder->writeBlockHeader(&stats, &header_length);
  assert(header_length == 4); /* 2 bytes for sentinel and amount*/
  byte* block = &data[0];
  encoder->encodeRange(block, block + amount);
  encoder->EndContextBlock();
  encoder->Finish();
  encoder->out_->write48bits(0xFF, len_pos);
  
  int bytes;
  uint64 packed_int = utils::packInteger(0x0F, &bytes);
  assert(bytes == 1);
  encoder->writePackedInteger(packed_int);
  delete encoder;

  bwtc::Decoder decoder(test_fname);
  char preproc = decoder.readGlobalHeader();
  assert(preproc == 'n');
  std::vector<uint64> context_lengths;
  uint64 compr_len = decoder.readBlockHeader(&context_lengths);
  assert(compr_len == 0xFF);
  assert(context_lengths.size() == 1);
  decoder.source_->Start();
  for(unsigned i = 0; i < data.size(); ++i) {
    assert(data[i] == decoder.decodeByte());
  }
  decoder.EndContextBlock();
  uint64 read_int = decoder.readPackedInteger();
  assert(read_int == packed_int);
  assert(bwtc::UnpackInteger(read_int) == 0x0F);
  assert(decoder.in_->CompressedDataEnding());
}
#endif

} //namespace tests

int main(int argc, char **argv) {
  if(argc < 2) {
    std::cout << "Fail: Give test file as a first argument.\n";
    return 1;
  }
  tests::test_fname = argv[1];
  tests::TestArithmeticCoding('n');
  tests::TestWritingAndReadingPackedIntegerList();
  /*
  tests::TestPackingIntegers(1000);
  tests::TestBlockArithmeticCoding(1000, 'n');
  //tests::TestWritingAndReadingHeadersAndSimpleData();
  //tests::TestWritingAndReadingBlockHeaders();
  */
  std::cout << "Encoder and Decoder passed all tests.\n";
}

