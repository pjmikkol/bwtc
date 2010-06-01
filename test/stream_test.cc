#include <cassert>

#include <string>
#include <vector>

#include <boost/filesystem/operations.hpp>
namespace fs = boost::filesystem;

#include "../stream.h"

namespace tests {

void BlockWriteTest();
void WriteToFile(const std::string& fname, int characters);
void TestFileSize(const std::string& fname, int characters);

const std::string fname("temptest.txt");


/* Test case for simple writing with the OutStream-object*/
void BlockWriteTest() {
  for (int i = 1; i < 10; i++) {
    int characters = 2000*i*i;
    WriteToFile(fname, characters);
    TestFileSize(fname, characters);
  }
}

void WriteToFile(const std::string& fname, int characters) {
  bwtc::OutStream f(fname);
  std::vector<char> data(characters,'a');
  f.WriteBlock(data.begin(), data.end());
}

void TestFileSize(const std::string& fname, int characters) {
  fs::path file(fname);
  assert(fs::file_size(file) == characters);
}

} //namespace tests


int main() {
  tests::BlockWriteTest();
  return 0;
}


