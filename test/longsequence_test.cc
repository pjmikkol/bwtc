/**
 * @file longsequence_test.cc
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
 * Testing for long sequences and border calculations.
 */

#include <cassert>
#include <cstdlib>
#include <ctime>

#include <iostream>
#include <string>
#include <vector>

#include "../preprocessors/test_preprocessor.h"
#include "../preprocessors/longsequences.h"
#include "../preprocessors/postprocessor.h"
#include "../globaldefs.h"
#include "../BlockManager.hpp"

#undef NDEBUG

namespace bwtc {
int verbosity = 3;
}

using bwtc::uint64;
using bwtc::byte;

namespace tests {

/******************** Testing calculation of borders ******************/
void ValidateBorder(int border, const std::vector<byte>& data) {
  for(int i = 0; i < border; ++i) {
    assert(data[i] == data[data.size() - border + i]);
  }
  for(int i = border + 1; static_cast<unsigned>(i) < data.size() - 1; ++i) {
    bool is_border = true;
    for(int j = 0; j < i; ++j) {
      if(data[j] != data[data.size() - i + j]) {
        is_border = false;
        break;
      }
    }
    assert(!is_border);
  }
}

void TestBorderWithRandom(int times, unsigned size_of_alphabet) {
  using bwtc::long_sequences::Border;
  if (size_of_alphabet > 256) size_of_alphabet = 256;
  srand(time(NULL));
  for(int j = 0; j < times; ++j) {
    int size_of_test = (rand() & 0xFFFF) + 2;
    std::vector<byte> data(size_of_test);
    for(int i = 0; i < size_of_test; ++i) {
      data[i] = rand() % size_of_alphabet;
    }
    Border<byte> b(data.size());
    int border = b(&data[0]);
    ValidateBorder(border, data);
  }
}

void TestWithGivenString(const std::string& str, int border) {
  using bwtc::long_sequences::Border;
  Border<char> b(str.size());
  assert(b(str.c_str()) == border);
}

void TestBorder() {
  TestBorderWithRandom(5,256);
  TestBorderWithRandom(5,128);
  TestBorderWithRandom(5,4);
  TestBorderWithRandom(5,2);
  TestWithGivenString(std::string("ababa"), 3);
  TestWithGivenString(std::string("ababababababa"), 11);
  TestWithGivenString(std::string("abbBbbab"), 2);
  TestWithGivenString(std::string("daadaa"), 3);
  TestWithGivenString(std::string("abbabaaabba"), 4);
}

void TestSequenceCompression(const std::string& source_name, int times,
                             uint64 block_size, unsigned window_size,
                             int threshold, int mem_constr)
{
  bwtc::BlockManager bm(block_size, 1);
  bwtc::TestPreProcessor pp(block_size);
  pp.AddBlockManager(&bm);
  pp.Connect(source_name);
  pp.InitializeTarget();
  uint64 data_size = pp.FillBuffer();
  std::vector<byte> original(pp.curr_block_->m_filled);
  std::copy(pp.curr_block_->begin(), pp.curr_block_->end(), original.begin());
  uint64 compressed_size;
  for(int j = 0; j < times; ++j) {
    compressed_size = bwtc::CompressSequences(pp.curr_block_->begin(),
                                              pp.curr_block_->m_filled,
                                              mem_constr, window_size,
                                              threshold);
    pp.curr_block_->m_filled = compressed_size;
  }
  uint64 uncompressed_size;
  for(int j = 0; j < times; ++j) {
    uncompressed_size = bwtc::UncompressSequences(pp.curr_block_->m_block,
                                                  pp.curr_block_->m_filled);
    pp.curr_block_->m_filled = uncompressed_size;
  }
  std::vector<byte>& uncompressed = *pp.curr_block_->m_block;

  assert(uncompressed_size == original.size());
  for(uint64 j = 0; j < data_size; ++j) {
    assert(uncompressed[j] == original[j]);
  }
}

} //namespace tests

int main(int argc, char **argv) {
  using namespace tests;
  TestBorder();
  uint64 block_size = 209715200;
  int times = 1;
  unsigned window_size = 16;
  int threshold = 128;
  int mem_constr = 2;
  if(argc > 2) window_size = atoi(argv[2]);
  if (argc > 3) threshold = atoi(argv[3]);
  if (argc > 4) mem_constr = atoi(argv[4]);
  if(argc > 5) times =  atoi(argv[5]);
  if (argc > 6) block_size = atoi(argv[6]);
  if(argc == 1) return 0;
  TestSequenceCompression(std::string(argv[1]), times, block_size, window_size,
                          threshold, mem_constr);
}
