#include <cassert>

#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <algorithm>

#include "stream.h"

namespace bwtc {


OutStream::OutStream(std::string file_name) :
    name_(file_name), target_(NULL), outfile_(NULL),
    buffer_(kDefaultBufferSize)
{
  if (name_ != "") {
    outfile_ = new std::ofstream(name_.c_str());
    target_ = new std::ostream_iterator<char>(*outfile_);
  } else
    target_ = new std::ostream_iterator<char>(std::cout);
  assert(target_);
}


OutStream::~OutStream() {
  if (outfile_) outfile_->close();
  delete outfile_;
  delete target_;
}


void OutStream::WriteBlock(std::vector<char>::const_iterator begin,
                           std::vector<char>::const_iterator end) {
  std::copy(begin, end, *target_);
}


} //namespace bwtc
