/**
 * @file BWTBlock.cpp
 * @author Pekka Mikkola <pmikkol@gmail.com>
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
 * Implementation of BWT-block.
 */

#include "BWTBlock.hpp"
#include "globaldefs.hpp"
#include "Streams.hpp"

#include <vector>

namespace bwtc {

BWTBlock::BWTBlock()
    : m_begin(0), m_length(0), m_isTransformed(true) {}

BWTBlock::BWTBlock(byte *data, uint32 length, bool isTransformed)
    : m_begin(data), m_length(length), m_isTransformed(isTransformed) {}

BWTBlock::BWTBlock(const BWTBlock& b)
    : m_begin(b.m_begin), m_length(b.m_length),
      m_LFpowers(b.m_LFpowers), m_isTransformed(b.m_isTransformed) {}

BWTBlock& BWTBlock::operator=(const BWTBlock& b) {
  m_begin = b.m_begin;
  m_length = b.m_length;
  m_LFpowers = b.m_LFpowers;
  m_isTransformed = b.m_isTransformed;
  return *this;
}

void BWTBlock::setBegin(byte* begin) {
  m_begin = begin;
}

void BWTBlock::setSize(uint32 length) {
  m_length = length;
}

size_t BWTBlock::writeHeader(OutStream* out) const {
  if(verbosity > 2) {
    std::clog << "Writing " << m_LFpowers.size() << " starting points."
              << std::endl;
  }
  size_t bytes = 1;
  byte s = (byte)(m_LFpowers.size()-1);
  out->writeByte(s);
  int bitsLeft = 8;
  for(size_t i = 0; i < m_LFpowers.size(); ++i) {
    for(int j = 30; j >= 0; --j) {
      s = (s << 1) | ((m_LFpowers[i] >> j) & 0x1);
      --bitsLeft;
      if(bitsLeft == 0) {
        out->writeByte(s);
        bitsLeft = 8;
        ++bytes;
      }
    }
  }
  if(bitsLeft < 8) {
    out->writeByte(s << bitsLeft);
    ++bytes;
  }
  return bytes;
}

void BWTBlock::readHeader(InStream* in) {
  uint32 LFpows = in->readByte()+1;
  if(verbosity > 2) {
    std::clog << "Reading " << LFpows << " starting points."
              << std::endl;
  }
  m_LFpowers.resize(LFpows);
  for(uint32 i = 0; i < LFpows; ++i) {
    uint32 pos = 0;
    for(uint32 j = 0; j < 31; ++j)
      pos = (pos << 1) | (in->readBit() ? 1 : 0);
    m_LFpowers[i] = pos;
  }
  in->flushBuffer();
}

void BWTBlock::prepareLFpowers(uint32 startingPoints) {
  if(m_length <= 256 || startingPoints == 0) m_LFpowers.resize(1);
  else if(startingPoints <= 256) m_LFpowers.resize(startingPoints);
  else m_LFpowers.resize(256);
}

} //namespace bwtc
