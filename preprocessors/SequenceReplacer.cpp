/**
 * @file SequenceReplacer.cpp
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
 * Implementation for preprocessor which replaces the often occurring long
 * strings.
 */
#include "SequenceReplacer.hpp"
#include "../Utils.hpp"
#include "../globaldefs.hpp"

#include <algorithm>
#include <vector>

namespace bwtc {

uint64 SequenceReplacer::initHash(const byte* data) const {
  
}

void SequenceReplacer::analyseData(const byte *data, size_t length, bool reset) {
  if(reset) resetAnalyseData();
  assert(m_phase <= 1);
  if(m_phase == 0) {
    resizeAndInitTable(length/m_windowSize);
  } 
  m_phase = 1;
  size_t mask = m_hashValues.size()-1;
  // mask == 2^k - 1 for some k
  assert((mask & (mask + 1)) == 0);
  initHashConstant();
  uint64 h = initHash(data);
  
  
}

void SequenceReplacer::initHashConstant() {
  m_hashRemovalConstant = 1;
  for(size_t i = 1; i < m_windowSize; ++i) {
    m_hashRemovalConstant *= s_hashConstant;
  }
}

void SequenceReplacer::resizeAndInitTable(size_t preference) {
  size_t s = (preference > 0)?utils::mostSignificantBit(preference):2;
  if(s >= (1 << 31)) {
    m_hashValues.resize(1 << 31);
  } else {
    m_hashValues.resize(s);
  }
  std::fill(m_hashValues.begin(), m_hashValues.end(), std::make_pair(0,0));
}

} //namespace bwtc
