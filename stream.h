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
  virtual ~OutStream();
  /* Writes first amount_of_bits of char to stream */
  void Write(char bits, int amount_of_bits);
  /* Writes chars in range [begin, end) to stream */
  void WriteBlock(std::vector<char>::const_iterator begin,
                  std::vector<char>::const_iterator end);
  void Flush();

  
 private:
  std::string name_;
  std::ostream_iterator<char>* target_;
  std::ofstream* outfile_;
  std::vector<char> buffer_;
  
  OutStream& operator=(const OutStream& os);
  OutStream(const OutStream& os);
};


} //namespace bwtc

#endif
