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
      m_useEscaping(sr.m_useEscaping),m_hashValues(sr.m_hashValues),
      m_sequences(sr.m_sequences), m_buckets(sr.m_buckets)
{
  std::copy(sr.m_frequencies, sr.m_frequencies + 256, m_frequencies);
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
  const uint32 halfWindow = (m_windowSize+1)/2;

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
    // TODO: ensure that if !foundDuplicates -> maxPos is not duplicate
    bool foundDuplicates = false;
    bool acceptedDuplicates = false;
    bool foundProper = false;
    for(size_t i = 1; i < bufferIndex; ++i) {
      size_t hash = buffer[i].first & mask;
      if(!foundDuplicates && hash == prevHash) {
        uint32 periodLength = buffer[i].second - buffer[i-1].second;
        if(periodLength >= halfWindow) {
          uint32 extraHash = buffer[i].first >> logMask;
          if(m_hashValues[hash].second != extraHash) continue;
          m_sequences.push_back(std::make_pair(hash, buffer[i-1].second));
          m_sequences.push_back(std::make_pair(hash, buffer[i].second));
          m_hashValues[hash].first += 2;
          calculateFrequencies(data, pos, buffer[i].second + periodLength);
          pos = buffer[i].second + periodLength;
          foundDuplicates = true;
          acceptedDuplicates = true;
          break;
        } else {
          if(maxPos == i-1) {
            maxCount = 0;
          }
          foundDuplicates = true;
          acceptedDuplicates = false;
          prevHash = hash;
          continue;
        }
      }
      prevHash = hash;
      size_t count = m_hashValues[prevHash].first;
      if(count >= maxCount) {
        maxCount = count;
        maxPos = i;
        foundProper = true;
      }
    }
    if(!acceptedDuplicates && foundProper) {
      size_t hash = buffer[maxPos].first;
      std::pair<uint32, uint32>& pair = m_hashValues[hash & mask];
      m_sequences.push_back(std::make_pair(hash & mask, buffer[maxPos].second));
      ++pair.first;
      pair.second = hash >> logMask;
      calculateFrequencies(data, pos, buffer[maxPos].second + m_windowSize);
      pos = buffer[maxPos].second + m_windowSize;
    } else if (!acceptedDuplicates && !foundProper) {
      pos += limit + 1;
    }
  }
}

void SequenceReplacer::sortIntoBuckets() {
  uint32 total = 0;
  for(size_t i = 0; i < m_hashValues.size(); ++i) {
    uint32 t = m_hashValues[i].first;
    if(t > 1) {
      m_hashValues[i].first = total; total += t;
    } else {
      m_hashValues[i].first = s_errorVal;
    }
  }
  if(m_verbose) {
    std::clog << "Total recurring sequences detected: " << total << std::endl;
  }

  m_buckets.resize(total);
  const uint32 size = m_sequences.size();
  for(uint32 i = 0; i < size; ++i) {
    std::pair<uint32, uint32> p = m_sequences[i];
    m_sequences[i].second |= 0x80000000;
    if(m_hashValues[p.first].first == s_errorVal) {
      continue;
    }
    p.second = i;
    m_buckets[m_hashValues[p.first].first++] = p;
  }
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
bool comparePosition(const std::pair<uint32, uint32>& p1,
                     const std::pair<uint32, uint32>& p2)
{
  return p1.second < p2.second;
}
} //namespace long_sequences

void SequenceReplacer::sortPositions(int begin, int end) {
  std::sort(&m_buckets[begin], &m_buckets[end],
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
    if(len == 1) m_buckets[begin].second |= 0x80000000;
    return;
  }
  /*if(len <= s_insertionSortLimit) {
    if(len > 1) insertionSort(begin, end, data);
    m_buckets[begin].second |= 0x80000000;
    return;
    }*/
  int pivot = end-1;
  int i = begin - 1, j = begin;
  for(; j < end - 1; ++j) {
    if(strCmp(m_sequences[m_buckets[pivot].second].second  & 0x7fffffff,
              m_sequences[m_buckets[j].second].second & 0x7fffffff, data) > 0)
    {
      ++i;
      std::swap(m_buckets[i], m_buckets[j]);
    }
  }
  sortSubBucket(begin, i+1, data);
  int equalBucket = j = i + 1;
  for(; j < end - 1; ++j) {
    if(strCmp(m_sequences[m_buckets[pivot].second].second & 0x7fffffff,
              m_sequences[m_buckets[j].second].second & 0x7fffffff, data) == 0)
    {
      ++i;
      std::swap(m_buckets[i], m_buckets[j]);
    }
  }
  ++i; std::swap(m_buckets[i], m_buckets[j]); ++i;
  sortPositions(equalBucket, i);
  m_buckets[equalBucket].second |= 0x80000000;
  sortSubBucket(i, end, data);
  m_buckets[begin].second |= 0x80000000;
}

void SequenceReplacer::sortAndMarkBuckets(const byte* data) {
  //TODO: fix if size == 0
  int prevPos = 0;
  size_t prevHash = m_buckets[0].first;
  int s = m_buckets.size();
  for(int i = 1; i < s; ++i) {
    if(prevHash != m_buckets[i].first) {
      sortSubBucket(prevPos, i, data);
      prevPos = i;
      prevHash = m_buckets[i].first;
    }
  }
  if(prevPos + 1 < s) {
    sortSubBucket(prevPos, s, data);
  }
  m_buckets[prevPos].second |= 0x80000000;
}

void SequenceReplacer::finishAnalysation() {
  assert(m_phase == 1);
}

void SequenceReplacer::nameRange(uint32 begin, uint32 end, uint32 name) {
  for(uint32 j = begin; j < end; ++j) {
    m_buckets[j].first = name;
  }
}

uint32 SequenceReplacer::nameHashValues() {
  uint32 name = 0;
  uint32 prev = 0;
  m_buckets[0].second &= 0x7fffffff;
  uint32 j = 0;
  for(uint32 i = 1; i < m_buckets.size(); ++i) {
    if(m_buckets[i].second & 0x80000000) {
      if(i > prev + 1) {
        nameRange(prev, i, name);
        ++name;
        for(uint32 k = prev; k < i; ++k) {
          m_buckets[j] = m_buckets[k];
          std::pair<uint32, uint32>& p = m_sequences[m_buckets[j].second];
          p.first = j++;
          p.second &= 0x7fffffff;
        }
      }
      prev = i;
      m_buckets[i].second &= 0x7fffffff;
    }
  }
  if(m_buckets.size() > prev + 1) {
    nameRange(prev, m_buckets.size(), name);
    for(uint32 k = prev; k < m_buckets.size(); ++k) {
      m_buckets[j] = m_buckets[k];
      std::pair<uint32, uint32>& p = m_sequences[m_buckets[j].second];
      p.first = j++;
      p.second &= 0x7fffffff;
    }
  }
  m_buckets.resize(j);
  return ++name;
}

bool SequenceReplacer::
validateRange(uint32 begin, uint32 end, const byte* data) const {
  if(m_buckets[begin].first == s_errorVal) {
    for(uint32 i = begin; i < end; ++i) {
      if(strCmp(m_sequences[m_buckets[i].second].second,
                m_sequences[m_buckets[i+1].second].second, data) == 0) {
        std::cout << "Discarded two equal strings." << std::endl;
        return false;
      }
    }
    return true;
  }
  if(begin + 1 == end) {
    if(m_buckets[begin].first != s_errorVal) {
      std::cout << "Wrong name for bucket of size 1" << std::endl;
      return false;
    }
    return true;
  }
  for(uint32 i = begin; i < end-1; ++i) {
    if(strCmp(m_sequences[m_buckets[i].second].second,
              m_sequences[m_buckets[i+1].second].second, data) != 0) {
      std::cout << "Same name with different strings" << std::endl;
      return false;
    }
    if(m_buckets[i].second >= m_buckets[i+1].second) {
      std::cout << "Positions " << m_buckets[i].second << " and "
                << m_buckets[i+1].second << " with same strings" << std::endl;
      return false;
    }
    if(m_sequences[m_buckets[i].second].first != i) {
      std::cout << "Incorrect pointer form sequences to buckets." << std::endl;
    }
  }
  return true;
}

bool SequenceReplacer::validatePhase2(const byte* data) const {
  uint32 begin = 0;
  for(uint32 i = 1; i < m_buckets.size(); ++i) {
    if(m_buckets[i].first != m_buckets[begin].first) {
      if(!validateRange(begin, i, data)) return false;
      if(strCmp(m_buckets[i].second, m_buckets[begin].second, data) == 0) {
        std::cout << "Different names with same strings" << std::endl;
        return false;
      }
      begin = i;
    }
  }
  return validateRange(begin, m_buckets.size(), data);
}

void SequenceReplacer::removeDeletedSequences() {
  size_t j = 0;
  for(size_t i = 0; i < m_sequences.size(); ++i) {
    if((m_sequences[i].second & 0x80000000) == 0) {
      m_sequences[j] = m_sequences[i];
      m_buckets[m_sequences[j].first].second = j;
      ++j;
    }
  }
  m_sequences.resize(j);
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
  m_hashValues.resize(separateStrings);
  if(m_verbose) {
    std::clog << separateStrings << " separate strings in " <<
        m_buckets.size() << " values." << std::endl;
  }

  if(validatePhase2(data)) {
    std::cout << "Everything ok" << std::endl;
  } else {
    std::cout << "Failure" << std::endl;
  }

  removeDeletedSequences();
  std::cout << m_sequences.size() << " " << m_buckets.size() << std::endl;
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
