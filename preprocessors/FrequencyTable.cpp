/**
 * @file FrequencyTable.cpp
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
 * Implementation of FrequencyTable-class which is used for bookkeeping
 * the frequencies of bytes.
 */

#include "../globaldefs.hpp"
#include "FrequencyTable.hpp"

#include <algorithm>
#include <cassert>
#include <utility>

namespace bwtc {

template <typename T>
bool comparePairSecondAsc(const T& p1, const T& p2) {
  return (p1.second < p2.second);
}

FrequencyTable::FrequencyTable() : m_last(256) {
  for(int i = 0; i < 256; ++i) {
    m_frequencies[i] = std::make_pair(static_cast<byte>(i), 0);
  }
  initLocations();
  assert(test());
}

FrequencyTable::FrequencyTable(const FrequencyTable& freqTable)
    : m_last(freqTable.m_last) {
  std::copy(freqTable.m_frequencies, freqTable.m_frequencies+256, m_frequencies);
  std::copy(freqTable.m_location, freqTable.m_location+256, m_location);
}

FrequencyTable& FrequencyTable::operator=(const FrequencyTable& freqTable) {
  m_last = freqTable.m_last;
  std::copy(freqTable.m_frequencies, freqTable.m_frequencies+256, m_frequencies);
  std::copy(freqTable.m_location, freqTable.m_location+256, m_location);
  return *this;
}


FrequencyTable::FrequencyTable(uint32* frequencies) : m_last(256) {
  initialize(frequencies);
}

void FrequencyTable::initialize(uint32* frequencies) {
  /* Assumes that frequencies has length of 256 */
  for(int i = 0; i < 256; ++i) {
    m_frequencies[i] = std::make_pair(static_cast<byte>(i), frequencies[i]);
  }
  std::sort(m_frequencies, m_frequencies + 256,
            comparePairSecondAsc<std::pair<byte, uint32> >);
  initLocations();
}


uint32 FrequencyTable::getFrequency(int i) const {
  assert(i >= 0 && i < m_last);
  return m_frequencies[i].second;
}

uint32 FrequencyTable::getFrequencyWithKey(byte key) const {
  return m_frequencies[m_location[key]].second;
}

byte FrequencyTable::getKey(int i) const {
  assert(i >= 0 && i <= 255);
  return m_frequencies[i].first;
}

bool FrequencyTable::decrease(byte key, uint32 value) {
  int freqIndex = m_location[key];
  assert(m_frequencies[freqIndex].second >= value);
  if(m_frequencies[freqIndex].second < value) return false;
  value = m_frequencies[freqIndex].second - value;
  std::pair<byte, uint32> pair =
      std::make_pair(m_frequencies[freqIndex].first, value);
  
  while (freqIndex > 0 && value < m_frequencies[freqIndex - 1].second)
  {
    ++m_location[m_frequencies[freqIndex - 1].first];
    m_frequencies[freqIndex] = m_frequencies[freqIndex - 1];
    --freqIndex;
  }
  m_frequencies[freqIndex] = pair;
  m_location[pair.first]= freqIndex;
  assert(m_frequencies[m_location[pair.first]].first == pair.first);
  return true;
}

void FrequencyTable::remove(byte key) {
  uint32 freqIndex = m_location[key];
  if(freqIndex < m_last) {
    --m_last;
    byte swapKey = m_frequencies[m_last].first;
    std::swap(m_frequencies[m_last], m_frequencies[freqIndex]);
    m_location[swapKey] = freqIndex;
    m_location[key] = m_last;
    increase(key, 0);
  }
}

void FrequencyTable::increase(byte key, uint32 value) {
  uint32 freqIndex = m_location[key];
  
  value += m_frequencies[freqIndex].second;
  assert(m_frequencies[freqIndex].first == key);
  std::pair<byte, uint32> pair = std::make_pair(key, value);
  
  while (freqIndex < m_last && value > m_frequencies[freqIndex + 1].second)
  {
    --m_location[m_frequencies[freqIndex + 1].first];
    m_frequencies[freqIndex] = m_frequencies[freqIndex + 1];
    ++freqIndex;
  }
  m_frequencies[freqIndex] = pair;
  m_location[key]= freqIndex;
  assert(m_frequencies[m_location[pair.first]].first == pair.first);
}

void FrequencyTable::initLocations() {
  for(int i = 0; i < 256; ++i) {
    m_location[m_frequencies[i].first] = i;
  }
}

bool FrequencyTable::test() {
  for(int i = 0; i < 256; ++i) {
    assert(m_frequencies[m_location[i]].first == i );
  }
  return true;
}


}
