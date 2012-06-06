/**
 * @file Streams.hpp
 * @author Pekka Mikkola <pjmikkol@cs.helsinki.fi>
 * @author Dominik Kempa <dominik.kempa@cs.helsinki.fi>
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
 * Header for InStream, OutStream, RawInStream and RawOutStream.
 */

#ifndef BWTC_STREAM_HPP_
#define BWTC_STREAM_HPP_

#include <cstdio>
#include <cassert>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

#include "globaldefs.hpp"

namespace bwtc {

/**
 * RawOutStream writes data into file or std::cout.
 *
 * On the compression pipeline RawOutStream-object is located at the end of
 * the pipeline. It is owned by Compressor-object which control the states of
 * stream. The actual writes are done by EntropyEncoder.
 *
 * RawOutStream is also at the end of the decompression pipeline right after
 * the postprocessors.
 *
 * @see Compressor
 */
class RawOutStream {
 public:
  explicit RawOutStream(const std::string& file_name);
  ~RawOutStream();

  /**
   * Writes given byte into target.
   *
   *@param b byte to be written
   */
  void writeByte(byte b);

  /**
   * Writes bytes from range [begin, end) to stream
   */
  void writeBlock(byte *begin, byte *end);

  bool isStdout() const { return m_fileptr == stdout; }

  long int getPos();
  void write48bits(uint64 to_written, long int position);
  void flush();

 private:
  static const uint32 kBufferSize = 1 << 16; // 64KB

  std::string m_name;
  FILE *m_fileptr;
  uint32 m_filled;
  byte *m_buffer;

  RawOutStream& operator=(const RawOutStream& os);
  RawOutStream(const RawOutStream& os);
};


struct RawInStream {
 public:
  explicit RawInStream(const std::string &file_name);
  ~RawInStream();

  /* Copies block from stream to given byte array.
   * Returns the number of read bytes. */
  uint64 readBlock(byte *to, uint64 max_block_size);

  inline bool readBit() {
    if (m_bitsInBuffer == 0) {
      m_buffer = static_cast<byte>(fetchByte());
      m_bitsInBuffer = 8;
    }
    return (m_buffer >> --m_bitsInBuffer) & 1;
  }

  inline byte readByte() {
    assert(m_bitsInBuffer < 8);
    byte nextByte = static_cast<byte>(fetchByte());
    m_buffer = (m_buffer << 8) | nextByte;
    return (m_buffer >> m_bitsInBuffer) & 0xff;
  }

  inline void flushBuffer() {
    m_bitsInBuffer = 0;
  }
  
  bool isStdin() const { return m_fileptr == stdin; }

  uint64 read48bits();

  bool compressedDataEnding() {
    /* Quick workaround. For some mysterious reason there is single
     * additional byte in the end of compressed file. It seems that
     * BitEncoder is responsible for this. */
    if (m_bigbuf_left > 1) return false;
    if (m_bigbuf_left == 0) {
      return (peekByte() == EOF);
    } else {
      return feof(m_fileptr);
    }
  }

 private:
  static const uint32 kBigbufSize = 1 << 16; // 64KB
 
  std::string m_name;
  FILE *m_fileptr;
  uint16 m_buffer;
  byte m_bitsInBuffer;
  byte *m_bigbuf;
  int32 m_bigbuf_left;
  uint32 m_bigbuf_pos;

  int fetchByte();
  int peekByte();
  RawInStream& operator=(const RawInStream& os);
  RawInStream(const RawInStream& os);
};

} //namespace bwtc

#endif
