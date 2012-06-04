/**
 * @file stream_test.cc
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
 * Testing of RawInStream and RawOutStream.
 *
 */

#include <cassert>

#include <iostream>
#include <string>
#include <vector>

#include <boost/filesystem/operations.hpp>
namespace fs = boost::filesystem;

#include "../globaldefs.hpp"
#include "../Streams.hpp"

#undef NDEBUG

using bwtc::uint64;
using bwtc::byte;

namespace tests {

std::string test_fname;

void BlockWriteTest();
void WriteToFile(const std::string& fname, unsigned characters);
void TestFileSize(const std::string& fname, unsigned characters);

void SimpleWriteReadTest() {
  bwtc::RawOutStream out(test_fname);
  out.writeByte('a');
  out.writeByte('b');
  out.writeByte('c');
  out.flush();
  bwtc::RawInStream in(test_fname);
  assert(in.readByte() == 'a');
  assert('b' == in.readByte());
  assert('c' == in.readByte());
}


/*********** begin: BlockFileWriteTest ***********/
void WriteToFile(const std::string& fname, unsigned characters) {
  bwtc::RawOutStream f(fname);
  std::vector<byte> data(characters,'a');
  f.writeBlock(&data[0], &data[0] + data.size());
}

void TestFileSize(const std::string& fname, unsigned characters) {
  fs::path file(fname);
  assert(fs::file_size(file) == characters);
}

/* Test case for simple writing with the RawOutStream-object*/
void BlockFileWriteTest() {
  for (int i = 1; i < 10; i++) {
    unsigned int characters = 2000*i*i;
    WriteToFile(test_fname, characters);
    TestFileSize(test_fname, characters);
  }
}
/*********** end: BlockFileWriteTest ***********/


void WriteToStreamTest() {
  bwtc::RawOutStream f("");
  std::string str = std::string("test\n");
  std::vector<byte> data(str.begin(), str.end());
  f.writeBlock(&data[0], &data[0] + data.size());
}

void EmptyWriteTest() {
  bwtc::RawOutStream f(test_fname);
  std::vector<byte> data(0);
  f.writeBlock(&data[0], &data[0] + data.size());
}

/*********** begin: ReadFromFileTest ***********/
void WriteAndRead(long fsize, long bsize) {  
  bwtc::RawOutStream o(test_fname);
  std::vector<byte> data(fsize, 'b');
  o.writeBlock(&data[0], &data[0] + data.size());
  o.flush();
  bwtc::RawInStream f(test_fname);
  std::vector<byte> block(bsize);
  long total = 0;
  while (std::streamsize read = f.readBlock(&block[0], 1000)) total += read;
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


int main(int argc, char **argv) {
  if(argc < 2) {
    std::cout << "Fail: Give test file as a first argument.\n";
    return 1;
  }
  tests::test_fname = argv[1];  
  tests::BlockFileWriteTest();
  tests::WriteToStreamTest();
  tests::EmptyWriteTest();
  tests::SimpleWriteReadTest();
  tests::ReadFromFileTest();
  std::cout << "Streams passed all tests.\n";
  return 0;
}
