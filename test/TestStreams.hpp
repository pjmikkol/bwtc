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

#include <vector>

namespace bwtc {
namespace tests {

class TestStream : public InStream, public OutStream {
 public:
  TestStream() : InStream(""), OutStream(""), m_currByte(0), m_currBit(0) {}
  
  void reset() { m_currBit = m_currByte = 0;  }

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
  
 private:
  std::vector<byte> m_data;
  uint32 m_currByte;
  uint32 m_currBit;
};

}} //namespace bwtc::tests




#endif
