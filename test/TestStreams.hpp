/**
 * @file TestStreams.hpp
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
 * Streams used in tests.
 */

#ifndef BWTC_TEST_STREAMS_HPP_
#define BWTC_TEST_STREAMS_HPP_

#include "../Streams.hpp"

#include <algorithm>
#include <vector>

namespace bwtc {
namespace tests {

class TestStream : public InStream, public OutStream {
 public:
  TestStream(std::vector<byte>& data)
      :  m_data(data), m_currByte(0), m_currBit(0) {}

  void reset() { m_currBit = m_currByte = 0;  }

  ~TestStream() {}
 
  bool readBit() {
    bool bit = (m_data[m_currByte] >> (7 - m_currBit++)) & 1;
    if(m_currBit == 8) {
      m_currBit = 0;
      ++m_currByte;
    }
    return bit;
  }

  void writeByte(byte b) {
    if(m_currByte == m_data.size()) {
      m_data.push_back(b);
      ++m_currByte;
    } else {
      m_data[m_currByte++] = b;
    }
  }

  void writeBlock(const byte* begin, const byte* end) {
    if(m_currBit) {
      ++m_currByte;
      m_currBit = 0;
    }
    m_data.resize((end - begin) + m_currByte);
    std::copy(begin, end, &m_data[m_currByte]);
  }

  long int getPos() {
    return m_currByte;
  }

  void flush() {}

  size_t readBlock(byte *to, size_t maxBlockSize) {
    size_t end = std::min(m_data.size(), m_currByte + maxBlockSize);
    std::copy(&m_data[m_currByte], &m_data[end], to);
    uint64 result = end - m_currByte;
    m_currByte = end;
    return result;
  }

  byte readByte() {
    return m_data[m_currByte++];
  }

  void flushBuffer() {
    if(m_currBit) {
      m_currBit = 0;
      ++m_currByte;
    }
  }

  void write48bits(uint64 toWritten, long int position) {
    position += 5;
    for(int i = 5; i >= 0; --i) {
      byte b = (toWritten >> i*8) & 0xff;
      m_data[position - i] = b;
    }
  }

  uint64 read48bits() {
    uint64 res = 0;
    for(size_t i = 0; i <= 5; ++i) {
      res <<= 8;
      res |= m_data[m_currByte++];
    }
    return res;
  }

  bool compressedDataEnding() {
    return m_currByte >= m_data.size();
  }
  
 private:
  std::vector<byte>& m_data;
  size_t m_currByte;
  uint32 m_currBit;
};

}} //namespace bwtc::tests




#endif
