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
#include <map>

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
      calculateFrequencies(data, pos, pos + limit + 1);
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

std::pair<uint32, uint32> SequenceReplacer::
findLeftLimit(uint32 sequence, uint32 offset) const {
  if(sequence == offset)
    return std::make_pair(m_sequences[sequence].second, offset);
  uint32 prev = sequence - 1 - offset;
  uint32 prevName = m_buckets[m_sequences[prev].first].first;
  uint32 prevLen = m_hashValues[prevName].second;
  if(prevLen > m_windowSize && (m_sequences[prev].second & 0x80000000) == 0) {
    return std::make_pair(m_sequences[sequence].second - m_sequences[prev].second - prevLen, offset);
  } else {
    return std::make_pair(m_sequences[sequence].second - m_sequences[prev].second, offset+1);
  }
}

std::pair<uint32, uint32> SequenceReplacer::
findRightLimit(uint32 sequence, uint32 offset, size_t length) const {
  if(sequence + offset + 1 == m_sequences.size())
    return std::make_pair(length - m_sequences[sequence].second - m_windowSize, offset);
  uint32 next = sequence + 1 + offset;
  uint32 nextName = m_buckets[m_sequences[next].first].first;
  uint32 nextLen = m_hashValues[nextName].second;
  if(nextLen > m_windowSize && (m_sequences[next].second & 0x80000000) == 0) {
    return std::make_pair(m_sequences[next].second - m_sequences[sequence].second - m_windowSize, offset);
  } else {
    return std::make_pair(m_sequences[next].second + m_windowSize - m_sequences[sequence].second, offset+1);
  }
}

uint32 SequenceReplacer::
expandToLeft(const byte* data, const std::vector<uint32>& elements,
             std::vector<std::pair<uint32, uint32> >& leftLimit) {
  uint32 expandToLeft = 0;
  bool expanding = true;
  while(expanding) {
    uint32 oPos = m_sequences[m_buckets[elements.back()].second].second;
    byte b = data[oPos - expandToLeft - 1];
    for(uint32 i = elements.size() - 1; i > 0; --i) {
      std::pair<uint32, uint32> leftL = leftLimit[i];
      if(leftL.first <= expandToLeft) {
        std::pair<uint32, uint32> nl = findLeftLimit(m_buckets[elements[i]].second, leftL.second);
        if(nl == leftL) {
          expanding = false;
          break;
        }
        leftLimit[i] = nl;
      }

      uint32 nPos = m_sequences[m_buckets[elements[i-1]].second].second;
      if(nPos < expandToLeft + 1 || nPos + m_windowSize + expandToLeft >= oPos) {
        expanding = false;
        break;
      }
      byte nb = data[nPos - expandToLeft - 1];
      if(nb != b) {
        expanding = false;
        break;
      }
      oPos = nPos;
    }
    std::pair<uint32, uint32> leftL = leftLimit[0];
    if(leftL.first <= expandToLeft) {
      std::pair<uint32, uint32> nl = findLeftLimit(m_buckets[elements[0]].second, leftL.second);
      if(nl == leftL) break;
      leftLimit[0] = nl;
    }
    if(expanding) ++expandToLeft;
  }
  return expandToLeft;
}

uint32 SequenceReplacer::
expandToRight(const byte* data, size_t length, const std::vector<uint32>& elements,
              std::vector<std::pair<uint32, uint32> >& rightLimit, uint32 leftExp) {
  uint32 expandToRight = 0;
  bool expanding = true;
  while(expanding) {
    uint32 oPos = m_sequences[m_buckets[elements[0]].second].second;
    byte b = data[oPos + m_windowSize + expandToRight];
    for(uint32 i = 0; i < elements.size() - 1; ++i) {
      std::pair<uint32, uint32> rightL = rightLimit[i];
      if(rightL.first <= expandToRight) {
        std::pair<uint32, uint32> rl = findRightLimit(m_buckets[elements[i]].second, rightL.second, length);
        if(rl == rightL) {
          expanding = false;
          break;
        }
        rightLimit[i] = rl;
      }

      uint32 nPos = m_sequences[m_buckets[elements[i+1]].second].second;
      if(nPos + m_windowSize + expandToRight >= length || oPos + m_windowSize + expandToRight + leftExp >= nPos) {
        expanding = false;
        break;
      }
      byte nb = data[nPos + m_windowSize + expandToRight];
      if(nb != b) {
        expanding = false;
        break;
      }
      oPos = nPos;
    }
    std::pair<uint32, uint32> rightL = rightLimit.back();
    if(rightL.first <= expandToRight) {
      std::pair<uint32, uint32> rl = findRightLimit(m_buckets[elements.back()].second, rightL.second, length);
      if(rl == rightL) break;
      rightLimit.back() = rl;
    }
    if(expanding) ++expandToRight;
  }
  return expandToRight;
}

void SequenceReplacer::
removeOverlappingSequences(const std::vector<uint32>& elements, uint32 leftOffset,
                           uint32 lengthOfSequence,
                           long_sequences::MaxHeap<uint32, uint32>& heap)
{
  std::map<uint32, uint32> updates;
  for(size_t i = 0; i < elements.size(); ++i) {
    int seq = m_buckets[elements[i]].second;
    m_sequences[seq].second -= leftOffset;
    uint32 pos = m_sequences[seq].second;
    for(int lseq = seq-1; lseq >= 0; --lseq) {
      if(m_sequences[lseq].second & 0x80000000) continue;
      std::pair<uint32, uint32>& p = m_sequences[lseq];
      uint32 llen = m_hashValues[m_buckets[p.first].first].second;
      if(p.second + llen > pos) {
        uint32 name = m_buckets[p.first].first;
        --m_hashValues[name].first;
        ++updates[name];
        m_buckets[p.first].second |= 0x80000000;
        p.second |= 0x80000000;
      } else {
        break;
      }
    }
    for(int rseq = seq+1; rseq < m_sequences.size(); ++rseq) {
      if(m_sequences[rseq].second & 0x80000000) continue;
      std::pair<uint32, uint32>& p = m_sequences[rseq];
      if(pos + lengthOfSequence > p.second) {
        uint32 name = m_buckets[p.first].first;
        --m_hashValues[name].first;
        ++updates[name];
        m_buckets[p.first].second |= 0x80000000;
        p.second |= 0x80000000;
      } else {
        break;
      }
    }
  }
  for(std::map<uint32, uint32>::const_iterator it = updates.begin();
      it != updates.end(); ++it) {
    heap.decrease(it->first, it->second);
  }
}

void SequenceReplacer::
expandStringsInBucket(uint32 name, uint32 begin, uint32 end,
                      const byte *data, size_t length,
                      long_sequences::MaxHeap<uint32, uint32>& heap) {
  std::vector<uint32> notDeleted;
  // First member of pair tells how many chars we can at least expand.
  // Second member tells offset to the previously considered sequence
  std::vector<std::pair<uint32, uint32> >leftLimit, rightLimit;
  for(uint32 i = begin; i < end; ++i) {
    if((m_buckets[i].second & 0x80000000) == 0) {
      notDeleted.push_back(i);
      std::pair<uint32,uint32> leftL = findLeftLimit(m_buckets[i].second, 0);
      if(leftL.first == 0) return;
      std::pair<uint32,uint32> rightL = findRightLimit(m_buckets[i].second, 0, length);
      if(rightL.first == 0) return;
      leftLimit.push_back(leftL);
      rightLimit.push_back(rightL);
    }
  }
  if(notDeleted.size() == 0) return;
  else if(notDeleted.size() == 1) {
    m_sequences[m_buckets[notDeleted[0]].second].second |= 0x80000000;
    m_buckets[notDeleted[0]].second |= 0x80000000;
    return;
  }

  uint32 expandedToLeft = expandToLeft(data, notDeleted, leftLimit);
  uint32 expandedToRight = expandToRight(data, length, notDeleted, rightLimit, expandedToLeft);
  uint32 len = expandedToLeft + expandedToRight + m_windowSize;

  if(len > m_windowSize) {
    removeOverlappingSequences(notDeleted, expandedToLeft, len, heap);
    m_hashValues[m_buckets[begin].first].second = len;
  }

  std::cout << "expanded " << expandedToLeft + expandedToRight << std::endl;
}

void SequenceReplacer::
validateCorrectOrder(uint32 begin, uint32 end) {
  for(uint32 i = begin; i < end-1; ++i) {
    bool fail = false;
    if((m_buckets[i].second & 0x7fffffff) >= (m_buckets[i+1].second & 0x7fffffff)) {
      std::cout << "Wrong order in buckets ";
      fail = true;
    }
    uint32 p1 = m_sequences[m_buckets[i].second & 0x7fffffff].second & 0x7fffffff;
    uint32 p2 = m_sequences[m_buckets[i+1].second & 0x7fffffff].second & 0x7fffffff;
    if(p1 >= p2) {
      std::cout << "Wrong order in sequences ";
      fail = true;
    }
    if(fail) {
      std::cout << std::endl;
      return;
    }
  }
}

void SequenceReplacer::decideLengths(const byte *data, size_t length) {
  long_sequences::MaxHeap<uint32, uint32> heap(m_hashValues.size());
  std::map<uint32, uint32> bucketLengths;
  for(uint32 i = 0; i < m_buckets.size(); ) {
    uint32 name = m_buckets[i].first;
    uint32 j = i+1;
    uint32 minLength = m_windowSize;
    while(j < m_buckets.size() && name == m_buckets[j].first) {
      uint32 len = m_sequences[m_buckets[j-1].second].second -
          m_sequences[m_buckets[j].second].second;
      if(len < minLength) minLength = len;
      ++j;
    }
    if(minLength < m_windowSize) m_hashValues[name].second = minLength;
    else {
      assert(j- i == m_hashValues[name].first);
      heap.insert(name, j - i, i);
      bucketLengths[name] = j - i;
    }
    i = j;
  }
  heap.prepare();

  while(!heap.empty()) {
    long_sequences::MaxHeap<uint32, uint32>::HeapElement el = heap.removeMax();
    std::map<uint32, uint32>::iterator it = bucketLengths.find(el.id);
    assert(it != bucketLengths.end());
    uint32 bucketLength = it->second;
    bucketLengths.erase(it);

    validateCorrectOrder(el.value, el.value + bucketLength);

    expandStringsInBucket(el.id, el.value, el.value + bucketLength, data, length, heap);
      
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
      std::cout << "Incorrect pointer from sequences to buckets." << std::endl;
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

// Removes deleted sequences from m_sequences and counts separate sequences
void SequenceReplacer::prepareForLengthAnalysation() {
  std::fill(m_hashValues.begin(), m_hashValues.end(),
            std::make_pair(0, m_windowSize));
  
  size_t j = 0;
  for(size_t i = 0; i < m_sequences.size(); ++i) {
    if((m_sequences[i].second & 0x80000000) == 0) {
      m_sequences[j] = m_sequences[i];
      uint32 posInBuckets = m_sequences[j].first;
      m_buckets[posInBuckets].second = j;
      ++m_hashValues[m_buckets[posInBuckets].first].first;
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

  prepareForLengthAnalysation();
  assert(m_sequences.size() ==  m_buckets.size());
  decideLengths(data, length);
  std::cout << "Decided lengths" << std::endl;
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
