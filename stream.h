#ifndef BWTC_STREAM_H_
#define BWTC_STREAM_H_

#include <iostream>
#include <iterator>
#include <string>
#include <vector>

namespace bwtc {

class OutStream {
 public:
  // TODO: find justified buffer size
  static const int kDefaultBufferSize = (1 << 12);
  explicit OutStream(std::string file_name);
  ~OutStream();
  /* Writes rightmost amount_of_bits of char to stream */
  void WriteBits(char bits, int amount_of_bits);
  /* Writes chars in range [begin, end) to stream */
  void WriteBlock(std::vector<char>::const_iterator begin,
                  std::vector<char>::const_iterator end);

 private:
  std::string name_;
  std::ostream_iterator<char>* to_;
  std::ofstream* outfile_;
  std::vector<char> buffer_;
  // TODO: work out WriteBits
  //std::vector<char>::iterator current_byte;
  //  int bits_left_in
  
  OutStream& operator=(const OutStream& os);
  OutStream(const OutStream& os);
};

class InStream {
 public:
  explicit InStream(std::string file_name);
  ~InStream();

 private:
  std::string name_;
  std::istream_iterator<char>* from_;
  std::ifstream* infile_;

  InStream& operator=(const InStream& os);
  InStream(const InStream& os);
};

} //namespace bwtc

#endif
