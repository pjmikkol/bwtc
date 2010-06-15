#include <cassert>

#include <iostream>
#include <string>
#include <vector>

#include <boost/filesystem/operations.hpp>
namespace fs = boost::filesystem;

#include "../globaldefs.h"
#include "../stream.h"
#include "testdefs.h"

namespace tests {

void BlockWriteTest();
void WriteToFile(const std::string& fname, unsigned characters);
void TestFileSize(const std::string& fname, unsigned characters);

void SimpleWriteReadTest() {
  bwtc::OutStream out(test_fname);
  out.WriteByte('a');
  out.WriteByte('b');
  out.WriteByte('c');
  out.Flush();
  bwtc::InStream in(test_fname);
  assert(in.ReadByte() == 'a');
  assert('b' == in.ReadByte());
  assert('c' == in.ReadByte());
}


/*********** begin: BlockFileWriteTest ***********/
void WriteToFile(const std::string& fname, unsigned characters) {
  bwtc::OutStream f(fname);
  std::vector<byte> data(characters,'a');
  f.WriteBlock(data.begin(), data.end());
}

void TestFileSize(const std::string& fname, unsigned characters) {
  fs::path file(fname);
  assert(fs::file_size(file) == characters);
}

/* Test case for simple writing with the OutStream-object*/
void BlockFileWriteTest() {
  for (int i = 1; i < 10; i++) {
    unsigned int characters = 2000*i*i;
    WriteToFile(test_fname, characters);
    TestFileSize(test_fname, characters);
  }
}
/*********** end: BlockFileWriteTest ***********/


void WriteToStreamTest() {
  bwtc::OutStream f("");
  std::string str = std::string("test\n");
  std::vector<byte> data(str.begin(), str.end());
  f.WriteBlock(data.begin(), data.end());  
}

void EmptyWriteTest() {
  bwtc::OutStream f(test_fname);
  std::vector<byte> data(0);
  f.WriteBlock(data.begin(), data.end());  
}

/*********** begin: ReadFromFileTest ***********/
void WriteAndRead(long fsize, long bsize) {  
  bwtc::OutStream o(test_fname);
  std::vector<byte> data(fsize, 'b');
  o.WriteBlock(data.begin(), data.end());
  o.Flush();
  bwtc::InStream f(test_fname);
  std::vector<byte> block(bsize);
  long total = 0;
  while (std::streamsize read = f.ReadBlock(&block[0], 1000)) total += read;
  assert(total == fsize);
}

void ReadFromFileTest() {
  for(int i = 0; i < 10; ++i) {
    for(int j = 1; j <= 10; ++j) {
      /* Numbers here are "random" */
      WriteAndRead(i*9677L, j*j*3244L);
    }
  }
}
/*********** end: ReadFromFileTest ***********/


} //namespace tests


int main() {
  tests::BlockFileWriteTest();
  tests::WriteToStreamTest();
  tests::EmptyWriteTest();
  tests::SimpleWriteReadTest();
  std::cout << "OutStream passed all tests.\n";
  tests::ReadFromFileTest();
  std::cout << "InStream passed all tests.\n";
  return 0;
}


