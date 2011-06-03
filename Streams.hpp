/**
 * @file Streams.hpp
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
 * Header for InStream and OutStream.
 */

#ifndef BWTC_STREAM_HPP_
#define BWTC_STREAM_HPP_

#include <iostream>
#include <iterator>
#include <string>
#include <vector>

#include "globaldefs.hpp"

namespace bwtc {

/**
 * OutStream writes data into file or std::cout.
 *
 * On the compression pipeline OutStream-object is located at the end of
 * the pipeline. It is owned by Encoder-object which feeds it with compressed
 * data.
 *
 * OutStream is also at the end of the decompression pipeline right after the
 * postprocessors.
 *
 * @see Encoder
 * @see PostProcessor
 */
class OutStream {
 public:
  explicit OutStream(std::string file_name);
  ~OutStream();

  /**
   * Writes given byte into target.
   *
   *@param b byte to be written
   */
  void WriteByte(byte b);

  /**
   * Writes bytes from range [begin, end) to stream
   *
   * @param begin iterator to the start of the range to be written
   * @param end iterator to the one past last of the range to be written
   */
  void WriteBlock(std::vector<byte>::const_iterator begin,
                  std::vector<byte>::const_iterator end);
  std::streampos GetPos() const; 
  void Write48bits(uint64 to_written, std::streampos position);
  void Flush();

 private:
  std::string m_name;
  std::ostream* m_to;
  std::ofstream* m_outfile;
  
  OutStream& operator=(const OutStream& os);
  OutStream(const OutStream& os);
};

class InStream {
 public:
  explicit InStream(std::string file_name);
  ~InStream();
  /* Copies block from stream to given char array.
   * Returns the number of read chars. */
  std::streamsize ReadBlock(byte* to, std::streamsize max_block_size);
  byte ReadByte() { return static_cast<byte>(m_from->get()); }
  uint64 Read48bits();
  bool CompressedDataEnding() {
    /* Quick workaround. For some mysterious reason there is single
     * additional byte in the end of compressed file. It seems that
     * BitEncoder is responsible for this. */
    if (m_from->eof()) return true;
    ReadByte();
    if (m_from->eof()) return true;
    m_from->unget();
    return false;
  }

 private:
  std::string m_name;
  std::istream* m_from;
  std::ifstream* m_infile;

  InStream& operator=(const InStream& os);
  InStream(const InStream& os);
};

} //namespace bwtc

#endif
