#include <cassert>

#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <algorithm>

#include "globaldefs.h"
#include "stream.h"

namespace bwtc {

OutStream::OutStream(std::string file_name) :
    name_(file_name), to_(NULL), outfile_(NULL),
    buffer_(kDefaultBufferSize)
{
  if (name_ != "") {
    outfile_ = new std::ofstream(name_.c_str());
    to_ = dynamic_cast<std::ostream*>(outfile_);
  } else
    to_ = &std::cout;
  assert(to_);
}

OutStream::~OutStream() {
  if (outfile_) {
    outfile_->close();
    delete outfile_;
  }
}

void OutStream::Flush() {
  to_->flush();
}

void OutStream::WriteByte(byte b) {
  to_->put(b);
}

std::streampos OutStream::GetPos() const {
  return to_->tellp();
}

void OutStream::WriteBlock(std::vector<byte>::const_iterator begin,
                           std::vector<byte>::const_iterator end) {
  std::copy(begin, end, std::ostream_iterator<byte>(*to_));
}

void OutStream::Write48bits(uint64 to_written, std::streampos position) {
  assert((to_written & (((uint64)0xFFFF) << 48)) == 0 );
  std::streampos current = to_->tellp();
  to_->seekp(position);
  for(int i = 5; i >= 0; --i) {
    byte b = 0xFF & (to_written >> i*8);
    to_->put(b);
  }
  to_->flush();
  to_->seekp(current);
}

uint64 InStream::Read48bits() {
  uint64 result = 0;
  for(int i = 0; i < 6; ++i) {
    result <<= 8;
    result |= ReadByte();
  }
  return result;
}

InStream::InStream(std::string file_name) :
    name_(file_name), from_(NULL), infile_(NULL)
{
  if (name_ != "") {
    infile_ = new std::ifstream(name_.c_str());
    from_ = dynamic_cast<std::istream*>(infile_);
  } else
    from_ = &std::cin;
  assert(from_);
}

InStream::~InStream() {
  if (infile_) {
    infile_->close();
    delete infile_;
  }
}

std::streamsize InStream::ReadBlock(byte* to, std::streamsize max_block_size) {
  if(! *from_) return 0; /* stream has reached EOF */
  from_->read(reinterpret_cast<char*>(to), max_block_size);
  return from_->gcount();
}

} //namespace bwtc
