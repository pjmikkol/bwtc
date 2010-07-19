#include "../preprocessors/longsequences.cc"

#include <cassert>
#include <cstdlib>
#include <ctime>

#include <iostream>
#include <string>
#include <vector>


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
  using long_sequences::Border;
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
  using long_sequences::Border;
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

} //namespace tests

int main() {
  using namespace tests;
  TestBorder();
}
