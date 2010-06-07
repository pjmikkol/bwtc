#include <cassert>

#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <algorithm>

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

//void OutStream::WriteBits(char bits, int amount_of_bits) {
  // if buffer_ is full WriteBlock(buffer_.begin(), buffer_.end())
//}

void OutStream::WriteBlock(std::vector<char>::const_iterator begin,
                           std::vector<char>::const_iterator end) {
  std::copy(begin, end, std::ostream_iterator<char>(*to_));
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

std::streamsize InStream::ReadBlock(std::vector<char>::iterator to,
                                    std::streamsize max_block_size) {
  if(! *from_) return 0;
  from_->read(&*to, max_block_size);
  return from_->gcount();
}

} //namespace bwtc
