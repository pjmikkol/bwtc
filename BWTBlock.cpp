/**
 * @file BWTBlock.cpp
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
 * Implementation of BWT-block.
 */

#include "BWTBlock.hpp"
#include "globaldefs.hpp"
#include "Streams.hpp"

#include <vector>

namespace bwtc {

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
}

void BWTBlock::writeHeader(OutStream* out) const {
  if(verbosity > 2) {
    std::clog << "Writing " << m_LFpowers.size() << " starting points."
              << std::endl;
  }
  //TODO: write them
}

void BWTBlock::prepareLFpowers(uint32 startingPoints) {
  if(m_length <= 256 || startingPoints == 0) m_LFpowers.resize(1);
  else if(startingPoints <= 256) m_LFpowers.resize(startingPoints);
  else m_LFpowers.resize(256);
}

} //namespace bwtc
