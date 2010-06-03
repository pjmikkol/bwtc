#include <cassert>

#include <iostream>
#include <string>
#include <vector>

#include <boost/filesystem/operations.hpp>
namespace fs = boost::filesystem;

#include "../stream.h"

namespace tests {

void BlockWriteTest();
void WriteToFile(const std::string& fname, unsigned characters);
void TestFileSize(const std::string& fname, unsigned characters);

const std::string fname("temptest.txt");

void WriteToFile(const std::string& fname, unsigned int characters) {
  bwtc::OutStream f(fname);
  std::vector<char> data(characters,'a');
  f.WriteBlock(data.begin(), data.end());
}

void TestFileSize(const std::string& fname, unsigned int characters) {
  fs::path file(fname);
  assert(fs::file_size(file) == characters);
}

/* Test case for simple writing with the OutStream-object*/
void BlockFileWriteTest() {
  for (int i = 1; i < 10; i++) {
    unsigned int characters = 2000*i*i;
    WriteToFile(fname, characters);
    TestFileSize(fname, characters);
  }
}

void WriteToStream() {
  bwtc::OutStream f("");
  std::string str = std::string("test\n");
  std::vector<char> data(str.begin(), str.end());
  f.WriteBlock(data.begin(), data.end());  
}

void EmptyWrite() {
  bwtc::OutStream f(fname);
  std::vector<char> data(0);
  f.WriteBlock(data.begin(), data.end());  
}

} //namespace tests


int main() {
  tests::BlockFileWriteTest();
  tests::WriteToStream();
  tests::EmptyWrite();
  std::cout << "OutStream passed all tests.\n";
  return 0;
}


