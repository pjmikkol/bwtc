#include <cassert>
#include <cstdlib>
#include <ctime>

#include <iostream>
#include <string>
#include <vector>

#include "../preprocessors/test_preprocessor.h"
#include "../preprocessors/longsequences.h"
#include "../globaldefs.h"
#include "../block_manager.h"

namespace bwtc {
int verbosity = 2;
}

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
    int border = Border(&data[0], data.size());
    ValidateBorder(border, data);
  }
}

void TestWithGivenString(const std::string& str, int border) {
  using bwtc::long_sequences::Border;
  assert(Border(str.c_str(), str.size()) == border);
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
                             int threshold)
{
  bwtc::BlockManager bm(block_size, 1);
  bwtc::TestPreProcessor pp(block_size);
  pp.AddBlockManager(&bm);
  pp.Connect(source_name);
  pp.InitializeTarget();
  uint64 data_size = pp.FillBuffer();
  //for(int j = 0; j < times; ++j) {
  bwtc::CompressSequences(pp.curr_block_->begin(), pp.curr_block_->filled_, 5,
                          window_size, threshold);
  bwtc::CompressSequences(pp.curr_block_->begin(), pp.curr_block_->filled_, 2,
                          window_size, threshold);
    //}
}

} //namespace tests

int main(int argc, char **argv) {
  using namespace tests;
  TestBorder();
  uint64 block_size = 209715200;
  int times = 1;
  unsigned window_size = 16;
  int threshold = 128;
  if(argc > 2) window_size = atoi(argv[2]);
  if (argc > 3) threshold = atoi(argv[3]);
  if(argc > 4) times =  atoi(argv[4]);
  if (argc > 5) block_size = atoi(argv[5]);
  if(argc == 1) return 0;
  TestSequenceCompression(std::string(argv[1]), times, block_size, window_size,
                          threshold);
}
