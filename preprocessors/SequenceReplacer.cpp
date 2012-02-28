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
#include "FrequencyTable.hpp"
#include "../Utils.hpp"
#include "../globaldefs.hpp"
#include "../Profiling.hpp"

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
    : m_hashValues(sr.m_hashValues), m_sequences(sr.m_sequences),
      m_buckets(sr.m_buckets), m_samples(sr.m_samples),
      m_hashRemovalConstant(sr.m_hashRemovalConstant), 
      m_dataLength(sr.m_dataLength), m_windowSize(sr.m_windowSize),
      m_numOfReplacements(sr.m_numOfReplacements),
      m_escapeByte(sr.m_escapeByte), m_phase(sr.m_phase),
      m_data(sr.m_data), m_verbose(sr.m_verbose),
      m_useEscaping(sr.m_useEscaping)
{
  std::copy(sr.m_frequencies, sr.m_frequencies + 256, m_frequencies);
}

SequenceReplacer::~SequenceReplacer() {}

void SequenceReplacer::resetAnalyseData() {
  std::fill(m_frequencies, m_frequencies + 256, 0);
  m_sequences.clear();
}

void SequenceReplacer::
deleteRemovedAndTakeSamples(std::vector<long_sequences::Sequence>& replaceables) {
  PROFILE("SequenceReplacer::deleteRemovedAndTakeSamples");
  using long_sequences::Sequence;
  uint32 i = 0;
  while(m_buckets[i].positionInPos & 0x80000000) ++i;

  uint32 name = m_buckets[i].name;
  std::pair<uint32, uint32>& p = m_hashValues[name];
  replaceables.push_back(Sequence(p.first, p.second, name, m_buckets[i].position));
  uint32 j = 0;
  m_buckets[j] = m_buckets[i];
  m_sequences[m_buckets[i].positionInPos].first = m_buckets[j].name;
  //m_sequences[m_buckets[i].positionInPos].first = j++;

  bool storedSample = true;
  for(++i; i < m_buckets.size(); ++i) {
    if(name != m_buckets[i].name) {
      storedSample = false;
      name = m_buckets[i].name;
    }
    if(!storedSample && (m_buckets[i].positionInPos & 0x80000000) == 0) {
      p = m_hashValues[name];
      //posInSeq = m_buckets[i].second;
      //replaceables.push_back(
      //    Sequence(p.first, p.second, name, m_sequences[posInSeq].second));
      replaceables.push_back(
          Sequence(p.first, p.second, name, m_buckets[i].position));
      storedSample = true;
    }
    if((m_buckets[i].positionInPos & 0x80000000) == 0) {
      m_buckets[j] = m_buckets[i];
      m_sequences[m_buckets[i].positionInPos].first = name;
      //m_sequences[m_buckets[i].second].first = j++;
    }
  }
  m_buckets.resize(j);  

  for(i = 0,j = 0; i < m_sequences.size(); ++i) {
    if((m_sequences[i].second & 0x80000000) == 0) {
      m_sequences[j] = m_sequences[i];
      //m_sequences[j].first = m_buckets[m_sequences[i].first].name;
      ++j;
    }
  }
  m_sequences.resize(j);
  assert(m_sequences.size() == m_buckets.size());
  m_buckets.clear();
  if(m_verbose) {
    std::clog << replaceables.size() << " different sequences to choose from."
              << std::endl;
  }
}

uint32 SequenceReplacer::
findReplaceableSequences(const std::vector<long_sequences::Sequence>& replaceables,
                         FrequencyTable& freqTable, uint32 maxReps) const
{
  PROFILE("SequenceReplacer::findReplaceableSequences");
  size_t current = 0;
  uint32 lim = std::min(maxReps, (uint32)replaceables.size());
  while(current < lim) {
    uint32 freqs[256] = {0};
    uint32 pos = replaceables[current].samplePosition;
    calculateFrequencies(pos, pos + replaceables[current].length, freqs);
    for(uint32 i = 0; i < 256; ++i) {
      if(freqs[i]) freqTable.decrease((byte)i, freqs[i]);
    }

    if(freqTable.getFrequency(current) + 1003 >=
       replaceables[current].count*(replaceables[current].length - 1)) {
      for(uint32 i = 0; i < 256; ++i) {
        if(freqs[i]) freqTable.increase((byte)i, freqs[i]);
      }
      break;
    }
    ++current;
  }
  return current;
}

uint32 SequenceReplacer::
findEscapeIndex(FrequencyTable& freqTable, uint32 freeSymbols,
                std::vector<long_sequences::Sequence>& replaceables,
                uint32 candidates) const
{
  if(candidates <= freeSymbols) return freeSymbols;
  int64 utility = 0;
  uint32 i = freeSymbols;
  for(;i < candidates; ++i) {
    utility += (replaceables[i].count*(replaceables[i].length-1) -
                freqTable.getFrequency(i) - replaceables[i].length - 5);
  }
  while(utility <= static_cast<int64>(freqTable.getFrequency(i)) && i > freeSymbols)
  {
    --i;
    uint32 freqs[256] = {0};
    uint32 pos = replaceables[i].samplePosition;
    calculateFrequencies(pos, pos + replaceables[i].length, freqs);
    for(uint32 j = 0; j < 256; ++j) {
      if(freqs[j]) freqTable.increase((byte)j, freqs[j]);
    }
    utility -= (replaceables[i].count*(replaceables[i].length-1) -
                freqTable.getFrequency(i) - replaceables[i].length - 5);
  }
  return i;
}

uint32 SequenceReplacer::decideReplacements() {
  if(m_buckets.size() == 0) return 0;
  assert(m_phase == 2);
  std::vector<long_sequences::Sequence> replaceables;
  deleteRemovedAndTakeSamples(replaceables);
  //TODO: try partial_sort
  std::sort(replaceables.rbegin(), replaceables.rend());

  FrequencyTable freqTable(m_frequencies);

  uint32 freeSymbols = 0;
  while(freqTable.getFrequency(freeSymbols) == 0) ++freeSymbols;

  uint32 escapeIndex;

  uint32 candidates = 0;
  if(m_useEscaping) {
    candidates = findReplaceableSequences(replaceables, freqTable, 254);
    escapeIndex = (candidates > freeSymbols)?
        findEscapeIndex(freqTable, freeSymbols, replaceables, candidates):freeSymbols;
  } else {
    if(freeSymbols > 0)
      candidates = findReplaceableSequences(replaceables, freqTable, freeSymbols);
    escapeIndex = freeSymbols;
  }

  m_numOfReplacements = (escapeIndex > freeSymbols)?
      escapeIndex: std::min(freeSymbols, candidates);
  if(m_verbose) {
    std::clog << "Replacing " << m_numOfReplacements << " sequences. ";
    if(m_numOfReplacements > freeSymbols)
      std::clog << "Made " << (escapeIndex - freeSymbols + 1)
                << " symbols free." << std::endl;
    else
      std::clog << "No symbols made free." << std::endl;
  }

  std::fill(m_isEscaped, m_isEscaped + 256, false);

  if(m_numOfReplacements > freeSymbols) {
    m_useEscaping = true;
    m_escapeByte = freqTable.getKey(escapeIndex);
    for(uint32 i = freeSymbols; i <= escapeIndex; ++i) {
      m_isEscaped[freqTable.getKey(i)] = true;
    }
  } else {
    m_useEscaping = false;
  }
  
  std::fill(m_hashValues.begin(), m_hashValues.end(), std::make_pair(0,0));

  for(uint i = 0; i < m_numOfReplacements; ++i) {
    uint32 name = replaceables[i].name;
    m_hashValues[name].first = freqTable.getKey(i);
    m_hashValues[name].second = replaceables[i].length;
    m_samples[name] = m_data + replaceables[i].samplePosition;
  }

  // To make writing phase a bit easier
  m_sequences.push_back(std::make_pair(0, m_dataLength));

  return m_numOfReplacements;
}

uint32 SequenceReplacer::writeAndPackInteger(byte* to, uint32 length) const {
  int packedLength;
  size_t packedInteger = utils::packInteger(length, &packedLength);
  for(int i = 0; i < packedLength; ++i) {
    byte toWritten = static_cast<byte>(packedInteger & 0xFF);
    to[i] = toWritten;
    packedInteger >>= 8;
  }
  return packedLength;
}

size_t SequenceReplacer::writeHeader(byte *to) const {
  if(m_numOfReplacements == 0) {
    to[0] = to[1] = 0;
  }
  uint32 j = 0; byte prev = 0;
  for(std::map<uint32, const byte*>::const_iterator it = m_samples.begin();
      it != m_samples.end(); ++it) {
    uint32 name = it->first;
    byte replacement = m_hashValues[name].first;
    to[j++] = replacement;
    prev = replacement;
    uint32 length = m_hashValues[name].second;
    j += writeAndPackInteger(to + j, length);
    std::copy(it->second, it->second + length, to + j);
    j += length;
  }
  to[j++] = prev;
  to[j++] = m_useEscaping?m_escapeByte:prev;
  return j;
}

uint32 SequenceReplacer::writeWithEscaping(const byte* begin, const byte* end,
                                           byte* dst) const {
  uint32 j = 0;
  for(; begin < end; ++begin) {
    if(m_isEscaped[*begin]) dst[j++] = m_escapeByte;
    dst[j++] = *begin;
  }
  return j;
}

size_t SequenceReplacer::
writeReplacedVersion(const byte *src, size_t length, byte *dst) const {
  PROFILE("SequenceReplacer::writeReplacedVersion");
  if(m_numOfReplacements == 0) {
    std::copy(src, src + length, dst);
    return length;
  }
  //Moved to decideReplacements
  //m_sequences.push_back(std::make_pair(0, length));
    
  uint32 currentSeq = 0;
  uint32 j = writeWithEscaping(src, src + m_sequences[currentSeq].second, dst);
  uint32 i = m_sequences[currentSeq].second;

  while(i < length) {
    uint32 name = m_sequences[currentSeq].first;
    if(m_hashValues[name].second > 0) {
      dst[j++] = m_hashValues[name].first;
      i += m_hashValues[name].second;
    }
    ++currentSeq;
    uint32 pos = m_sequences[currentSeq].second;
    j += writeWithEscaping(src + i, src + pos, dst + j);
    i = m_sequences[currentSeq].second;
  }
  return j;
}

uint64 SequenceReplacer::initHash(const byte* data) const {
  uint64 h = 0;
  for(size_t i = 0; i < m_windowSize; ++i) {
    h = h*s_hashConstant + data[i];
  }
  return h;
}

void SequenceReplacer::
calculateFrequencies(uint32 begin, uint32 end, uint32 *f) const {
  while(begin < end) ++f[m_data[begin++]];
}

struct buf_struct {
  buf_struct() : hash(0), extraHash(0), count(0) {}
  buf_struct(size_t h, uint32 eh, uint32 c)
      : hash(h), extraHash(eh), count(c) {}
  size_t hash;
  uint32 extraHash;
  uint32 count;
};

void SequenceReplacer::scanAndStore() {
  PROFILE("SequenceReplacer::scanAndStore");
  size_t mask = m_hashValues.size()-1;
  size_t logMask = utils::logFloor(mask+1);
  // mask == 2^k - 1 for some k
  assert((mask & (mask + 1)) == 0);

  uint32 pos = 0;
  const uint32 limit = 2 * m_windowSize - 1;
  const uint32 halfWindow = (m_windowSize+1)/2;

  std::vector<buf_struct> buffer;
  buffer.resize(m_windowSize);

  while(pos + limit < m_dataLength) {
    uint64 h = initHash(m_data + pos);

    //Store the pairs <hash, count> into the buffer
    //Two separate versions. The second one seems to work better on
    //non-monotonous inputs
#if 0
    buffer[0] = buf_struct(h, m_hashValues[h & mask].second, m_hashValues[h & mask].first);
    for(uint32 i = 1; i < m_windowSize; ++i) {
      byte next = m_data[pos + i - 1 + m_windowSize];
      h -= m_data[pos + i - 1]*m_hashRemovalConstant;
      h *= s_hashConstant;
      h += next;
      buffer[i] = buf_struct(h, m_hashValues[h & mask].second, m_hashValues[h & mask].first);
    }
#else    
    buffer[0] = buf_struct(h, 0, 0);
    for(uint32 i = 1; i < m_windowSize; ++i) {
      byte next = m_data[pos + i - 1 + m_windowSize];
      h -= m_data[pos + i - 1]*m_hashRemovalConstant;
      h *= s_hashConstant;
      h += next;
      buffer[i] = buf_struct(h, 0, 0);
    }
    for(uint32 i = 0; i < m_windowSize; ++i) {
      std::pair<uint32, uint32> p = m_hashValues[buffer[i].hash & mask];
      buffer[i].extraHash = p.second;
      buffer[i].count = p.first;
    }
#endif

    
    uint32 maxCount = buffer[0].count;
    uint32 maxPos = 0;
    bool foundSuitable = ((buffer[0].hash >> logMask) & 0xffffffff) == buffer[0].extraHash;
    for(uint32 i = 1; i < m_windowSize; ++i) {
      if(buffer[i].count != 0 && ((buffer[i].hash >> logMask) & 0xffffffff) != buffer[i].extraHash) continue;
      if(!foundSuitable || buffer[i].count >= maxCount) {
        maxCount = buffer[i].count;
        maxPos = i;
        foundSuitable = true;
      }
    }
    
    if(foundSuitable) {
      size_t hash = buffer[maxPos].hash;
      std::pair<uint32, uint32>& pair = m_hashValues[hash & mask];
      ++pair.first;
      pair.second = hash >> logMask;
      m_sequences.push_back(std::make_pair(hash & mask, pos + maxPos));
      calculateFrequencies(pos, pos + maxPos + m_windowSize, m_frequencies);
      pos += maxPos + m_windowSize;
    } else {
      calculateFrequencies(pos, pos + limit, m_frequencies);
      pos += limit;
    }
  }
}

void SequenceReplacer::sortIntoBuckets() {
  PROFILE("SequenceReplacer::sortIntoBuckets");
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
    m_buckets[m_hashValues[p.first].first++] = bucket_struct(p.first, p.second, i);
  }
}

int SequenceReplacer::
strCmp(uint32 pos1, uint32 pos2) const {
  for(uint32 i = 0; i < m_windowSize; ++i) {
    int diff = m_data[pos1++] - m_data[pos2++];
    if(diff != 0) return diff;
  }
  return 0;
}

namespace long_sequences {
bool comparePosition(const SequenceReplacer::bucket_struct& p1,
                     const SequenceReplacer::bucket_struct& p2)
{
  return p1.position < p2.position;
}
} //namespace long_sequences

void SequenceReplacer::sortPositions(int begin, int end) {
  std::sort(&m_buckets[begin], &m_buckets[end],
            long_sequences::comparePosition);
}

/* TODO: Recognize and mark subbuckets
void SequenceReplacer::insertionSort(int begin, int end, const byte* data) {
  assert(end - begin > 0);
  for(int i = begin + 1; i < end; ++i) {
    int j = i;
    std::pair<uint32, uint32> val = m_sequences[j];
    int cmp = 0;
    while(j > begin) {
      cmp = strCmp(m_sequences[j].second, m_sequences[j-1].second);
      if(cmp > 0 || (cmp == 0 && m_sequences[j].second > m_sequences[j-1].second))
        break;
      m_sequences[j] = m_sequences[j-1];
      --j;
    }
    m_sequences[j] = val;
  }
}
*/

void SequenceReplacer::
sortSubBucket(int begin, int end, bool positionsUnordered) {
  int len = end - begin;
  if(len <= 1) {
    if(len == 1) m_buckets[begin].positionInPos |= 0x80000000;
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
    if(strCmp(m_buckets[pivot].position, m_buckets[j].position) > 0) {
      ++i;
      std::swap(m_buckets[i], m_buckets[j]);
    }
  }
  sortSubBucket(begin, i+1, false);
  positionsUnordered = positionsUnordered || i + 1 != begin;
  j = i + 1;
  int equalBucket = j;
  for(; j < end - 1; ++j) {
    if(strCmp(m_buckets[pivot].position, m_buckets[j].position) == 0) {
      ++i;
      std::swap(m_buckets[i], m_buckets[j]);
    }
  }
  ++i; std::swap(m_buckets[i], m_buckets[j]); ++i;
  if(positionsUnordered) sortPositions(equalBucket, i);
  m_buckets[equalBucket].positionInPos |= 0x80000000;
  sortSubBucket(i, end, true);
  m_buckets[begin].positionInPos |= 0x80000000;
}

void SequenceReplacer::sortAndMarkBuckets() {
  PROFILE("SequenceReplacer::sortAndMarkBuckets");
  if(m_buckets.size() == 0) return;
  int prevPos = 0;
  size_t prevHash = m_buckets[0].name;
  int s = m_buckets.size();
  for(int i = 1; i < s; ++i) {
    if(prevHash != m_buckets[i].name) {
      sortSubBucket(prevPos, i, false);
      prevPos = i;
      prevHash = m_buckets[i].name;
    }
  }
  if(prevPos + 1 < s) {
    sortSubBucket(prevPos, s, false);
  }
  m_buckets[prevPos].positionInPos |= 0x80000000;
}

void SequenceReplacer::finishAnalysation() {
  assert(m_phase == 2);
}

void SequenceReplacer::nameRange(uint32 begin, uint32 end, uint32 name) {
  for(uint32 j = begin; j < end; ++j) { m_buckets[j].name = name; }
}

std::pair<uint32, uint32> SequenceReplacer::
findLeftLimit(uint32 sequence, uint32 offset) const {
  if(sequence == offset) {
    return std::make_pair(m_sequences[sequence].second, offset);
  }
  int prev = sequence - offset - 1;
  while(prev > 0 && (m_sequences[prev].second & 0x80000000)) {
    --prev; ++offset;
  }
  uint32 prevName = m_buckets[m_sequences[prev].first].name;
  uint32 prevLen = m_hashValues[prevName].second;
  if(prevLen > m_windowSize && (m_sequences[prev].second & 0x80000000) == 0) {
    return std::make_pair(m_sequences[sequence].second - m_sequences[prev].second - prevLen, offset);
  } else {
    return std::make_pair(m_sequences[sequence].second - m_sequences[prev].second, offset+1);
  }
}

std::pair<uint32, uint32> SequenceReplacer::
findRightLimit(uint32 sequence, uint32 offset) const {
  if(sequence + offset + 1 == m_sequences.size())
    return std::make_pair(m_dataLength - m_sequences[sequence].second - m_windowSize, offset);
  uint32 next = sequence + offset + 1;
  while(next < m_sequences.size() && (m_sequences[next].second & 0x80000000)) {
    ++offset; ++next;
  }
  uint32 nextName = m_buckets[m_sequences[next].first].name;
  uint32 nextLen = m_hashValues[nextName].second;
  if(nextLen > m_windowSize) {
    return std::make_pair(m_sequences[next].second - m_sequences[sequence].second - m_windowSize, offset-1);
  } else {
    return std::make_pair(m_sequences[next].second - m_sequences[sequence].second, offset);
  }
}

uint32 SequenceReplacer::
expandToLeft(const std::vector<uint32>& elements,
             std::vector<std::pair<uint32, uint32> >& leftLimit) {
  uint32 expandToLeft = 0;
  bool expanding = true;
  while(expanding) {
    uint32 oPos = m_buckets[elements.back()].position;
    byte b = m_data[oPos - expandToLeft - 1];
    for(uint32 i = elements.size() - 1; i > 0; --i) {
      std::pair<uint32, uint32> leftL = leftLimit[i];
      if(leftL.first <= expandToLeft) {
        std::pair<uint32, uint32> nl = findLeftLimit(m_buckets[elements[i]].positionInPos, leftL.second);
        if(nl.first <= expandToLeft) {
          expanding = false;
          break;
        }
        leftLimit[i] = nl;
      }

      uint32 nPos = m_buckets[elements[i-1]].position;
      if(nPos < expandToLeft + 1 || nPos + m_windowSize + expandToLeft >= oPos) {
        expanding = false;
        break;
      }
      byte nb = m_data[nPos - expandToLeft - 1];
      if(nb != b) {
        expanding = false;
        break;
      }
      oPos = nPos;
    }
    std::pair<uint32, uint32> leftL = leftLimit[0];
    if(leftL.first <= expandToLeft) {
      std::pair<uint32, uint32> nl = findLeftLimit(m_buckets[elements[0]].positionInPos, leftL.second);
      if(nl.first == expandToLeft) break;
      leftLimit[0] = nl;
    }
    if(expanding) ++expandToLeft;
  }
  return expandToLeft;
}

uint32 SequenceReplacer::
expandToRight(const std::vector<uint32>& elements,
              std::vector<std::pair<uint32, uint32> >& rightLimit, uint32 leftExp) {
  uint32 expandToRight = 0;
  bool expanding = true;
  while(expanding) {
    uint32 oPos = m_buckets[elements[0]].position;
    byte b = m_data[oPos + m_windowSize + expandToRight];
    for(uint32 i = 0; i < elements.size() - 1; ++i) {
      std::pair<uint32, uint32> rightL = rightLimit[i];
      if(rightL.first <= expandToRight) {
        std::pair<uint32, uint32> rl = findRightLimit(m_buckets[elements[i]].positionInPos, rightL.second);
        if(rl.first <= expandToRight) {
          expanding = false;
          break;
        }
        rightLimit[i] = rl;
      }

      uint32 nPos = m_buckets[elements[i+1]].position;
      if(nPos + m_windowSize + expandToRight >= m_dataLength || oPos + m_windowSize + expandToRight + leftExp >= nPos) {
        expanding = false;
        break;
      }
      byte nb = m_data[nPos + m_windowSize + expandToRight];
      if(nb != b) {
        expanding = false;
        break;
      }
      oPos = nPos;
    }
    std::pair<uint32, uint32> rightL = rightLimit.back();
    if(rightL.first <= expandToRight) {
      std::pair<uint32, uint32> rl = findRightLimit(m_buckets[elements.back()].positionInPos, rightL.second);
      if(rl.first <= expandToRight) break;
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
    int seq = m_buckets[elements[i]].positionInPos;
    m_buckets[elements[i]].position -= leftOffset;
    m_sequences[seq].second -= leftOffset;
    uint32 pos = m_sequences[seq].second;
    for(int lseq = seq-1; lseq >= 0; --lseq) {
      if((m_sequences[lseq].second & 0x80000000) != 0) continue;
      std::pair<uint32, uint32>& p = m_sequences[lseq];

      //uint32 llen = m_hashValues[m_buckets[p.first].name].second;
      //if(p.second + llen > pos) {
      if(p.second + m_windowSize > pos) {
        uint32 name = m_buckets[p.first].name;
        //--m_hashValues[name].first;
        ++updates[name];
        m_buckets[p.first].positionInPos |= 0x80000000;
        p.second |= 0x80000000;
      } else {
        break;
      }
    }
    for(int rseq = seq+1; rseq < (int)m_sequences.size(); ++rseq) {
      if((m_sequences[rseq].second & 0x80000000) != 0) continue;
      std::pair<uint32, uint32>& p = m_sequences[rseq];
      if(pos + lengthOfSequence > p.second) {
        uint32 name = m_buckets[p.first].name;
        //--m_hashValues[name].first;
        ++updates[name];
        m_buckets[p.first].positionInPos |= 0x80000000;
        p.second |= 0x80000000;
      } else {
        break;
      }
    }
  }
  for(std::map<uint32, uint32>::const_iterator it = updates.begin();
      it != updates.end(); ++it) {
    heap.decrease(it->first, it->second);
    m_hashValues[it->first].first -= it->second;
  }
}

void SequenceReplacer::
expandStringsInBucket(uint32 begin, uint32 end,
                      long_sequences::MaxHeap<uint32, uint32>& heap) {
  std::vector<uint32> notDeleted;
  // First member of pair tells how many chars we can at least expand.
  // Second member tells offset to the previously considered sequence
  std::vector<std::pair<uint32, uint32> >leftLimit, rightLimit;
  for(uint32 i = begin; i < end; ++i) {
    if((m_buckets[i].positionInPos & 0x80000000) == 0) {
      notDeleted.push_back(i);
      std::pair<uint32,uint32> leftL = findLeftLimit(m_buckets[i].positionInPos, 0);
      if(leftL.first == 0) return;
      std::pair<uint32,uint32> rightL = findRightLimit(m_buckets[i].positionInPos, 0);
      if(rightL.first == 0) return;
      leftLimit.push_back(leftL);
      rightLimit.push_back(rightL);
    }
  }
  if(notDeleted.size() == 0) return;
  else if(notDeleted.size() == 1) {
    m_sequences[m_buckets[notDeleted[0]].positionInPos].second |= 0x80000000;
    m_buckets[notDeleted[0]].positionInPos |= 0x80000000;
    return;
  }

  uint32 expandedToLeft = expandToLeft(notDeleted, leftLimit);
  uint32 expandedToRight = expandToRight(notDeleted, rightLimit, expandedToLeft);
  uint32 len = expandedToLeft + expandedToRight + m_windowSize;

  if(len > m_windowSize) {
    removeOverlappingSequences(notDeleted, expandedToLeft, len, heap);
    m_hashValues[m_buckets[begin].name].second = len;
  }
}

void SequenceReplacer::
validateCorrectOrder(uint32 begin, uint32 end) {
  for(uint32 i = begin; i < end-1; ++i) {
    bool fail = false;
    if((m_buckets[i].positionInPos & 0x7fffffff) >= (m_buckets[i+1].positionInPos & 0x7fffffff)) {
      std::cout << "Wrong order in buckets ";
      fail = true;
    }
    uint32 p1 = m_sequences[m_buckets[i].positionInPos & 0x7fffffff].second & 0x7fffffff;
    uint32 p2 = m_sequences[m_buckets[i+1].positionInPos & 0x7fffffff].second & 0x7fffffff;
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

void SequenceReplacer::decideLengths() {
  PROFILE("SequenceReplacer::decideLengths");
  long_sequences::MaxHeap<uint32, uint32> heap(m_hashValues.size());
  std::map<uint32, uint32> bucketLengths;
  for(uint32 i = 0; i < m_buckets.size(); ) {
    uint32 name = m_buckets[i].name;
    uint32 j = i+1;
    while(j < m_buckets.size() && name == m_buckets[j].name) {
      ++j;
    }
    assert(j- i == m_hashValues[name].first);
    heap.insert(name, j - i, i);
    bucketLengths[name] = j - i;
    i = j;
  }
  heap.prepare();

  while(!heap.empty()) {
    long_sequences::MaxHeap<uint32, uint32>::HeapElement el = heap.removeMax();
    std::map<uint32, uint32>::iterator it = bucketLengths.find(el.id);
    assert(it != bucketLengths.end());
    uint32 bucketLength = it->second;
    bucketLengths.erase(it);
    expandStringsInBucket(el.value, el.value + bucketLength, heap);
  }
}

uint32 SequenceReplacer::nameHashValues() {
  if(m_buckets.size() == 0) return 0;
  uint32 name = 0;
  uint32 prev = 0;
  m_buckets[0].positionInPos &= 0x7fffffff;
  uint32 j = 0;
  for(uint32 i = 1; i < m_buckets.size(); ++i) {
    if(m_buckets[i].positionInPos & 0x80000000) {
      if(i > prev + 1) {
        nameRange(prev, i, name);
        ++name;
        for(uint32 k = prev; k < i; ++k) {
          m_buckets[j] = m_buckets[k];
          std::pair<uint32, uint32>& p = m_sequences[m_buckets[j].positionInPos];
          p.first = j++;
          p.second &= 0x7fffffff;
        }
      }
      prev = i;
      m_buckets[i].positionInPos &= 0x7fffffff;
    }
  }
  if(m_buckets.size() > prev + 1) {
    nameRange(prev, m_buckets.size(), name);
    for(uint32 k = prev; k < m_buckets.size(); ++k) {
      m_buckets[j] = m_buckets[k];
      std::pair<uint32, uint32>& p = m_sequences[m_buckets[j].positionInPos];
      p.first = j++;
      p.second &= 0x7fffffff;
    }
  }
  m_buckets.resize(j);
  return ++name;
}

bool SequenceReplacer::
validateRange(uint32 begin, uint32 end) const {
  if(m_buckets[begin].name == s_errorVal) {
    for(uint32 i = begin; i < end; ++i) {
      if(strCmp(m_buckets[i].position, m_buckets[i+1].position) == 0) {
        std::cout << "Discarded two equal strings." << std::endl;
        return false;
      }
    }
    return true;
  }
  if(begin + 1 == end) {
    if(m_buckets[begin].name != s_errorVal) {
      std::cout << "Wrong name for bucket of size 1" << std::endl;
      return false;
    }
    return true;
  }
  for(uint32 i = begin; i < end-1; ++i) {
    if(strCmp(m_buckets[i].position, m_buckets[i+1].position) != 0 &&
       (m_buckets[i+1].positionInPos & 0x80000000) == 0 &&
       (m_buckets[i].positionInPos & 0x80000000) == 0) {
      std::cout << "Same name with different strings in indices " << i
                << " and " << i+1 << " in positions " << m_buckets[i].position
                << " and " << m_buckets[i+1].position << std::endl;
      for(uint32 k = m_buckets[i].position, j = 0; j < m_windowSize; ++j) {
        std::cout << m_data[k+j];
      }
      std::cout << " with name " << m_buckets[i].name << std::endl;
      for(uint32 k = m_buckets[i+1].position, j = 0; j < m_windowSize; ++j) {
        std::cout << m_data[k+j];
      }
      std::cout << " with name " << m_buckets[i+1].name << std::endl;
      return false;
    }
    if(m_buckets[i].position >= m_buckets[i+1].position &&
       (m_buckets[i+1].positionInPos & 0x80000000) == 0 &&
       (m_buckets[i].positionInPos & 0x80000000) == 0) {
      std::cout << "Positions " << m_buckets[i].position << " and "
                << m_buckets[i+1].position << " with same strings" << std::endl;
      return false;
    }
    if((m_buckets[i].positionInPos & 0x80000000) == 0 && m_sequences[m_buckets[i].positionInPos].first != i) {
      std::cout << "Incorrect pointer from sequences to buckets." << std::endl;
    }
  }
  return true;
}

bool SequenceReplacer::validatePhase2() const {
  uint32 begin = 0;
  for(uint32 i = 1; i < m_buckets.size(); ++i) {
    if(m_buckets[i].name != m_buckets[begin].name) {
      if(!validateRange(begin, i)) return false;
      if(strCmp(m_buckets[i].position, m_buckets[begin].position) == 0 &&
         (m_buckets[begin].positionInPos & 0x80000000) == 0 &&
         (m_buckets[i].positionInPos & 0x80000000) == 0) {
        std::cout << "Different names with same strings in indices "
                  << begin << " and " << i << " with names "
                  << m_buckets[begin].name << " and " << m_buckets[i].name
                  << std::endl;
        std::cout << "Lengths: " << m_hashValues[m_buckets[begin].name].second
                  << " and "
                  << m_hashValues[m_buckets[i].name].second << std::endl;
        for(uint32 k = m_buckets[i].position, j = 0; j < m_windowSize; ++j) {
          std::cout << m_data[k+j];
        }
        std::cout << std::endl;
        return false;
      }
      begin = i;
    }
  }
  return validateRange(begin, m_buckets.size());
}

// Removes deleted sequences from m_sequences and counts separate sequences
void SequenceReplacer::prepareForLengthAnalysation() {
  PROFILE("SequenceReplacer::prepareForLengthAnalysation");
  std::fill(m_hashValues.begin(), m_hashValues.end(), std::make_pair(0, m_windowSize));
  
  size_t j = 0;
  for(size_t i = 0; i < m_sequences.size(); ++i) {
    if((m_sequences[i].second & 0x80000000) == 0) {
      m_sequences[j] = m_sequences[i];

      uint32 posInBuckets = m_sequences[j].first;
      m_buckets[posInBuckets].positionInPos = j;
      ++m_hashValues[m_buckets[posInBuckets].name].first;
      ++j;
    }
  }
  m_sequences.resize(j);
}


void SequenceReplacer::analyseData(const byte *data, size_t length, bool reset) {
  m_data = data;
  m_dataLength = length;
  if(reset) resetAnalyseData();
  assert(m_phase <= 1);
  if(m_phase == 0) {
    resizeAndInitTable(length/m_windowSize);
  } 
  initHashConstant();

  assert(m_windowSize <= length);
  scanAndStore();

  sortIntoBuckets();
  sortAndMarkBuckets();

  uint32 separateStrings = nameHashValues();
  m_hashValues.resize(separateStrings);
  if(m_verbose) {
    std::clog << separateStrings << " separate strings in " <<
        m_buckets.size() << " values." << std::endl;
  }

  prepareForLengthAnalysation();
  assert(m_sequences.size() ==  m_buckets.size());
  decideLengths();
  m_phase = 2;
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
