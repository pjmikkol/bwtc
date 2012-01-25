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

SequenceReplacer::SequenceReplacer(bool useEscaping)
    : m_windowSize(s_defaultWindowSize), m_numOfReplacements(0),
      m_escapeByte(0), m_phase(0), m_verbose(false), m_useEscaping(useEscaping)
{}

SequenceReplacer::SequenceReplacer(bool useEscaping, bool verbose)
    : m_windowSize(s_defaultWindowSize), m_numOfReplacements(0),
      m_escapeByte(0), m_phase(0), m_verbose(verbose), m_useEscaping(useEscaping)
{}

SequenceReplacer::SequenceReplacer(const SequenceReplacer& sr)
    : m_windowSize(sr.m_windowSize), m_numOfReplacements(sr.m_numOfReplacements),
      m_escapeByte(sr.m_escapeByte), m_phase(sr.m_phase), m_verbose(sr.m_verbose),
      m_useEscaping(sr.m_useEscaping)

{
  std::copy(sr.m_frequencies, sr.m_frequencies + 256, m_frequencies);
  m_hashValues = sr.m_hashValues;
  m_sequences = sr.m_sequences;
}

SequenceReplacer::~SequenceReplacer() {}

void SequenceReplacer::resetAnalyseData() {
  std::fill(m_frequencies, m_frequencies + 256, 0);
  m_sequences.clear();
}

size_t SequenceReplacer::decideReplacements() {

}

size_t SequenceReplacer::writeHeader(byte *to) const {

}

size_t SequenceReplacer::
writeReplacedVersion(const byte *src, size_t length, byte *dst) const {

}

uint64 SequenceReplacer::initHash(const byte* data) const {
  uint64 h = 0;
  for(size_t i = 0; i < m_windowSize; ++i) {
    h = h*s_hashConstant + data[i];
  }
  return h;
}

void SequenceReplacer::
calculateFrequencies(const byte* data, uint32 begin, uint32 end)  {
  while(begin < end) ++m_frequencies[data[begin++]];
}

void SequenceReplacer::scanAndStore(const byte* data, size_t length) {
  size_t mask = m_hashValues.size()-1;
  size_t logMask = utils::logFloor(mask+1);
  // mask == 2^k - 1 for some k
  assert((mask & (mask + 1)) == 0);

  uint32 pos = 0;
  const uint32 limit = 2 * m_windowSize - 2;
  std::vector<std::pair<size_t, uint32> > buffer;
  buffer.resize(m_windowSize-1);
  //TODO: Is the use of buffer too slow?
  while(pos + limit < length) {
    uint32 bufferIndex = 0;
    uint64 h = initHash(data + pos);
    if(((h >> logMask)&0xffffffff) == m_hashValues[h & mask].second || m_hashValues[h & mask].first == 0 )
      buffer[bufferIndex++] = std::make_pair(h, pos);

    for(size_t i = 1; i < m_windowSize - 1; ++i) {
      h = (h - data[pos + i - 1]*m_hashRemovalConstant)*s_hashConstant + data[pos + i - 1 + m_windowSize];
      if(((h >> logMask)&0xffffffff) == m_hashValues[h & mask].second  || m_hashValues[h & mask].first == 0 )
        buffer[bufferIndex++] = std::make_pair(h, pos+i);
    }
    if(bufferIndex == 0) { pos += limit; continue; }

    std::sort(buffer.begin(), buffer.begin()+bufferIndex);

    size_t prevHash = buffer[0].first & mask;
    size_t maxCount = m_hashValues[prevHash].first;
    size_t maxPos = 0;
    bool duplicates = false;
    for(size_t i = 1; i < bufferIndex; ++i) {
      size_t hash = buffer[i].first & mask;
      if(hash == prevHash) {
        size_t last = gatherDuplicates(i-1, buffer, bufferIndex, hash, buffer[i].first >> logMask, mask);
        calculateFrequencies(data, pos, last + m_windowSize);
        pos = last + m_windowSize;
        duplicates = true;
        break;
      }
      prevHash = hash;
      size_t count = m_hashValues[prevHash].first;
      if(count > maxCount) {
        maxCount = count;
        maxPos = i;
      }
    }
    if(!duplicates) {
      size_t hash = buffer[maxPos].first;
      std::pair<uint32, uint32>& pair = m_hashValues[hash & mask];
      m_sequences.push_back(std::make_pair(hash & mask, buffer[maxPos].second));
      ++pair.first;
      pair.second = hash >> logMask;
      calculateFrequencies(data, pos, buffer[maxPos].second + m_windowSize);
      pos = buffer[maxPos].second + m_windowSize;
    }
  }
}

void SequenceReplacer::sortIntoBuckets() {
  uint32 total = 0;
  for(size_t i = 0; i < m_hashValues.size(); ++i) {
    uint32 t = m_hashValues[i].first;
    if(t > 1) {
      total += t; m_hashValues[i].first = total;
    } else {
      m_hashValues[i].first = s_errorVal;
    }
  }
  if(m_verbose) {
    std::clog << "Total recurring sequences detected: " << total << std::endl;
  }
  std::vector<std::pair<uint32, uint32> > sequences;
  sequences.resize(total);
  const size_t size = m_sequences.size();
  for(size_t i = 0; i < size; ++i) {
    std::pair<uint32, uint32> p = m_sequences.back();
    m_sequences.pop_back();
    if(m_hashValues[p.first].first == s_errorVal) continue;
    sequences[--m_hashValues[p.first].first] = p;
  }
  assert(m_sequences.size() == 0);
  std::swap(m_sequences, sequences);
}

int SequenceReplacer::
strCmp(uint32 pos1, uint32 pos2, const byte* data) const {
  for(uint32 i = 0; i < m_windowSize; ++i) {
    int diff = data[pos1++] - data[pos2++];
    if(diff != 0) return diff;
  }
  return 0;
}

namespace long_sequences {
int comparePosition(const std::pair<uint32, uint32>& p1,
                    const std::pair<uint32, uint32>& p2)
{
  return p1.second < p2.second;
}
} //namespace long_sequences

void SequenceReplacer::sortPositions(int begin, int end) {
  std::sort(&m_sequences[begin], &m_sequences[end],
            long_sequences::comparePosition);
}

// TODO: Recognize and mark subbuckets
void SequenceReplacer::insertionSort(int begin, int end, const byte* data) {
  assert(end - begin > 0);
  for(int i = begin + 1; i < end; ++i) {
    int j = i;
    std::pair<uint32, uint32> val = m_sequences[j];
    int cmp = 0;
    while(j > begin) {
      cmp = strCmp(m_sequences[j].second, m_sequences[j-1].second, data);
      if(cmp > 0 || (cmp == 0 && m_sequences[j].second > m_sequences[j-1].second))
        break;
      m_sequences[j] = m_sequences[j-1];
      --j;
    }
    m_sequences[j] = val;
  }
}


void SequenceReplacer::
sortSubBucket(int begin, int end, const byte* data) {
  int len = end - begin;
  if(len <= 1) {
    if(len == 1) m_sequences[begin].second |= 0x80000000;
    return;
  }
  /*if(len <= s_insertionSortLimit) {
    if(len > 1) insertionSort(begin, end, data);
    m_sequences[begin].second |= 0x80000000;
    return;
    }*/
  int pivot = end-1;
  int i = begin - 1, j = begin;
  for(; j < end - 1; ++j) {
    if(strCmp(m_sequences[pivot].second, m_sequences[j].second, data) > 0) {
      ++i;
      std::swap(m_sequences[i], m_sequences[j]);
    }
  }
  sortSubBucket(begin, i+1, data);
  int equalBucket = j = i + 1;
  for(; j < end - 1; ++j) {
    if(strCmp(m_sequences[pivot].second, m_sequences[j].second, data) == 0) {
      ++i;
      std::swap(m_sequences[i], m_sequences[j]);
    }
  }
  ++i; std::swap(m_sequences[i], m_sequences[j]); ++i;
  sortPositions(equalBucket, i);
  m_sequences[equalBucket].second |= 0x80000000;
  sortSubBucket(i, end, data);
  m_sequences[begin].second |= 0x80000000;
}

void SequenceReplacer::sortAndMarkBuckets(const byte* data) {
  //TODO: fix if size == 0
  int prevPos = 0;
  size_t prevHash = m_sequences[0].first;
  int s = m_sequences.size();
  for(int i = 1; i < s; ++i) {
    if(prevHash != m_sequences[i].first) {
      sortSubBucket(prevPos, i, data);
      prevPos = i;
      prevHash = m_sequences[i].first;
    }
  }
  if(prevPos + 1 < s) {
    sortSubBucket(prevPos, s, data);
  }
  m_sequences[prevPos].second |= 0x80000000;
}

void SequenceReplacer::finishAnalysation() {
  assert(m_phase == 1);
}

size_t SequenceReplacer::
gatherDuplicates(size_t index,
                 const std::vector<std::pair<size_t, uint32> >& buffer,
                 size_t bufferSize, size_t hash, size_t extraHash, size_t mask)
{
  std::pair<uint32, uint32>& pair = m_hashValues[hash];
  uint32 prevPos = buffer[index].second;
  ++pair.first;
  pair.second = extraHash;
  m_sequences.push_back(std::make_pair(hash, buffer[index++].second));
  for(; index < bufferSize; ++index) {
    if((buffer[index].first & mask) != hash) break;
    if(buffer[index].second < s_minPeriod + prevPos) continue;
    prevPos = buffer[index].second;
    ++pair.first;
    m_sequences.push_back(std::make_pair(hash, buffer[index].second));
  }
  return m_sequences.back().second;
}

void SequenceReplacer::nameRange(uint32 begin, uint32 end, uint32 name) {
  for(uint32 j = begin; j < end; ++j) {
    m_sequences[j].first = name;
  }
}

uint32 SequenceReplacer::nameHashValues() {
  uint32 name = 0;
  uint32 prev = 0;
  m_sequences[0].second &= 0x7fffffff;
  uint32 j = 0;
  for(uint32 i = 1; i < m_sequences.size(); ++i) {
    if(m_sequences[i].second & 0x80000000) {
      if(i > prev + 1) {
        nameRange(prev, i, name);
        ++name;
        for(uint32 k = prev; k < i; ++k) {
          m_sequences[j++] = m_sequences[k];
        }
      }
      prev = i;
      m_sequences[i].second &= 0x7fffffff;
    }
  }
  if(m_sequences.size() > prev + 1) {
    nameRange(prev, m_sequences.size(), name);
    for(uint32 k = prev; k < m_sequences.size(); ++k) {
      m_sequences[j++] = m_sequences[k];
    }
  }
  m_sequences.resize(j);
  return ++name;
}

bool SequenceReplacer::
validateRange(uint32 begin, uint32 end, const byte* data) const {
  if(m_sequences[begin].first == s_errorVal) {
    for(uint32 i = begin; i < end; ++i) {
      if(strCmp(m_sequences[i].second, m_sequences[i+1].second, data) == 0) {
        std::cout << "Discarded two equal strings." << std::endl;
        return false;
      }
    }
    return true;
  }
  if(begin + 1 == end) {
    if(m_sequences[begin].first != s_errorVal) {
      std::cout << "Wrong name for bucket of size 1" << std::endl;
      return false;
    }
    return true;
  }
  for(uint32 i = begin; i < end-1; ++i) {
    if(strCmp(m_sequences[i].second, m_sequences[i+1].second, data) != 0) {
      std::cout << "Same name with different strings" << std::endl;
      return false;
    }
    if(m_sequences[i].second >= m_sequences[i+1].second) {
      std::cout << "Positions " << m_sequences[i].second << " and "
                << m_sequences[i+1].second << " with same strings" << std::endl;
      return false;
    }
  }
  return true;
}

bool SequenceReplacer::validatePhase2(const byte* data) const {
  uint32 begin = 0;
  for(uint32 i = 1; i < m_sequences.size(); ++i) {
    if(m_sequences[i].first != m_sequences[begin].first) {
      if(!validateRange(begin, i, data)) return false;
      if(strCmp(m_sequences[i].second, m_sequences[begin].second, data) == 0) {
        std::cout << "Different names with same strings" << std::endl;
        return false;
      }
      begin = i;
    }
  }
  return validateRange(begin, m_sequences.size(), data);
}

void SequenceReplacer::analyseData(const byte *data, size_t length, bool reset) {
  if(reset) resetAnalyseData();
  assert(m_phase <= 1);
  if(m_phase == 0) {
    resizeAndInitTable(length/m_windowSize);
  } 
  std::cout << "Initted table" << std::endl;
  m_phase = 1;
  initHashConstant();
  std::cout << "Initted h-c" << std::endl;

  assert(m_windowSize <= length);
  scanAndStore(data, length);
  std::cout << "Scanned" << std::endl;

  m_phase = 2;
  sortIntoBuckets();
  std::cout << "Sorted into buckets" << std::endl;
  sortAndMarkBuckets(data);
  std::cout << "Sorted and marked buckets" << std::endl;
  uint32 separateStrings = nameHashValues();
  if(m_verbose) {
    std::clog << separateStrings << " separate strings in " <<
        m_sequences.size() << " values." << std::endl;
  }

  if(validatePhase2(data)) {
    std::cout << "Everything ok" << std::endl;
  } else {
    std::cout << "Failure" << std::endl;
  }
}

void SequenceReplacer::initHashConstant() {
  m_hashRemovalConstant = 1;
  for(size_t i = 1; i < m_windowSize; ++i) {
    m_hashRemovalConstant *= s_hashConstant;
  }
}

void SequenceReplacer::resizeAndInitTable(size_t preference) {
  size_t s = (preference > 0)?utils::mostSignificantBit(preference):2;
  if(s >= (1UL << 31)) {
    m_hashValues.resize(1UL << 31);
  } else {
    m_hashValues.resize(s);
  }
  std::fill(m_hashValues.begin(), m_hashValues.end(), std::make_pair(0,0));
}

} //namespace bwtc
