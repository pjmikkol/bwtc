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
    to_ = new std::ostream_iterator<char>(*outfile_);
  } else
    to_ = new std::ostream_iterator<char>(std::cout);
  assert(to_);
}


OutStream::~OutStream() {
  if (outfile_) outfile_->close();
  delete outfile_;
  delete to_;
}


//inline void OutStream::WriteBits(char bits, int amount_of_bits) {
  // if buffer_ is full WriteBlock(buffer_.begin(), buffer_.end())
//}

void OutStream::WriteBlock(std::vector<char>::const_iterator begin,
                           std::vector<char>::const_iterator end) {
  std::copy(begin, end, *to_);
}


InStream::InStream(std::string file_name) :
    name_(file_name), from_(NULL), infile_(NULL)
{
  if (name_ != "") {
    infile_ = new std::ifstream(name_.c_str());
    from_ = new std::istream_iterator<char>(*infile_);
  } else
    from_ = new std::istream_iterator<char>(std::cin);
  assert(from_);
}


InStream::~InStream() {
  if (infile_) infile_->close();
  delete infile_;
  delete from_;
}


} //namespace bwtc
