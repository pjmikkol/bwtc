/**
 * @file SequenceReplacer.hpp
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
 * Header for preprocessor which replaces the often occurring long strings.
 */
#ifndef SEQUENCE_REPLACER_HPP_
#define SEQUENCE_REPLACER_HPP_

#include "../globaldefs.hpp"

#include <cassert>
#include <vector>
#include <utility>

namespace bwtc {

namespace long_sequences {

template <typename K, typename V>
class MaxHeap {
 public:
  struct HeapElement {
    HeapElement(uint32 i, K k, V v) : id(i), key(k), value(v) {}
    uint32 id;
    K key;
    V value;
  };

  MaxHeap(uint32 maxIds) {
    m_positions.resize(maxIds);
    std::fill(m_positions.begin(),m_positions.end(), -1);
  }

  // Doesn't move elements. Need to call prepare after inserts
  void insert(size_t id, K key, V value) {
    int pos = m_positions[id];
    if(m_positions[id] >= 0) {
      m_heap[pos].key = key;
      m_heap[pos].value = value;
    } else {
      m_positions[id] = m_heap.size();
      m_heap.push_back(HeapElement(id, key, value));
    }
  }

  bool empty() {
    return m_heap.size() == 0;
  }
  
  HeapElement removeMax() {
    HeapElement el = m_heap[0];
    uint32 last = m_heap.size();
    if(last > 1) {
      --last;
      std::swap(m_heap[0], m_heap[last]);
      std::swap(m_positions[m_heap[0].id], m_positions[m_heap[last].id]);
    } 
    m_heap.pop_back();
    return el;
  }
  
  void decrease(uint32 id, K newKey) {
    m_heap[m_positions[id]].key = newKey;
    heapify(m_positions[id]);
  }
  
  void prepare() {
    initLocations();
    buildMaxHeap();
  }
  
 private:
  void initLocations() {
    for(uint32 i = 0; i < m_heap.size(); ++i) {
      m_positions[m_heap[i].id] = i;
    }
  }

#define parent(x) (((x)-1) >> 1)
#define left(x) (2*(x) + 1)
#define right(x) (2*(x) + 2)
  
  void heapify(int i) {
    int l = left(i), r = right(i);
    while(r < (int)m_heap.size()) {
      int largest = (m_heap[l].key < m_heap[r].key)? r : l;
      if(m_heap[i].key < m_heap[largest].key) {
        std::swap(m_positions[m_heap[largest].id], m_positions[m_heap[i].id]);
        std::swap(m_heap[largest], m_heap[i]);
        i = largest;
        l = left(i), r = right(i);
      } else {
        return;
      }
    }
    if(l+1 == (int)m_heap.size() && m_heap[i].key < m_heap[l].key) {
      std::swap(m_positions[m_heap[l].id], m_positions[m_heap[i].id]);
      std::swap(m_heap[l], m_heap[i]);
    }
  }

  void buildMaxHeap() {
    for(int i = parent(m_heap.size() - 1); i >= 0; --i) heapify(i);
  }

#undef parent
#undef left
#undef right
  std::vector<HeapElement> m_heap;
  // Mapping from id's to positions in heap
  std::vector<int> m_positions;  
};

} //namespace long_sequences

class SequenceReplacer {
 public:
  SequenceReplacer(bool useEscaping);
  SequenceReplacer(bool useEscaping, bool verbose);
  SequenceReplacer(const SequenceReplacer& sr);
  ~SequenceReplacer();

  void analyseData(const byte *data, size_t length, bool reset=true);

  void resetAnalyseData();

  void finishAnalysation();

  size_t decideReplacements();

  size_t writeHeader(byte *to) const;

  size_t writeReplacedVersion(const byte *src, size_t length, byte *dst) const;

 private:
  void resizeAndInitTable(size_t preference);
  uint64 initHash(const byte* data) const;
  void initHashConstant();
  void scanAndStore();
  void calculateFrequencies(const byte* data, uint32 begin, uint32 end);
  void sortIntoBuckets();
  void sortSubBucket(int begin, int end);
  int strCmp(uint32 pos1, uint32 pos2) const;
  void sortPositions(int begin, int end);
  void insertionSort(int begin, int end, const byte* data);
  void sortAndMarkBuckets();
  /**Names hash values and removes useless strings. Returns number of
   * separate string.*/
  uint32 nameHashValues();
  void nameRange(uint32 begin, uint32 end, uint32 name);
  bool validatePhase2() const;
  bool validateRange(uint32 begin, uint32 end) const;
  void prepareForLengthAnalysation();
  void decideLengths();
  void expandStringsInBucket(uint32 begin, uint32 end,
                             long_sequences::MaxHeap<uint32, uint32>& heap);
  std::pair<uint32, uint32> findLeftLimit(uint32 sequence, uint32 offset) const;
  std::pair<uint32, uint32> findRightLimit(uint32 sequence, uint32 offset) const;
  uint32 expandToLeft(const std::vector<uint32>& elements,
                      std::vector<std::pair<uint32, uint32> >& leftLimit);
  uint32 expandToRight(const std::vector<uint32>& elements,
                      std::vector<std::pair<uint32, uint32> >& rightLimit, uint32 leftExp);
  void removeOverlappingSequences(const std::vector<uint32>& elements, uint32 leftOffset,
                                  uint32 lengthOfSequence,
                                  long_sequences::MaxHeap<uint32, uint32>& heap);
  void validateCorrectOrder(uint32 begin, uint32 end);
  
  
  static const uint64 s_hashConstant = 37;
  static const uint32 s_errorVal = 0xffffffff;
  static const uint32 s_defaultWindowSize = 32;
  static const uint32 s_maxExpandableBucket = 10;
  //static const int s_insertionSortLimit = 10;

  /**Stores the frequencies of bytes. */
  size_t m_frequencies[256];

  /**Stores the hash values, counters and lengths of replaceable sequences.
   * When scanning, the first one is counter and second is h* */
  std::vector<std::pair<uint32, uint32> > m_hashValues;

  std::vector<std::pair<uint32, uint32> > m_sequences;

  std::vector<std::pair<uint32, uint32> > m_buckets;

  uint64 m_hashRemovalConstant;

  uint32 m_dataLength;

  uint32 m_windowSize;
  
  /** Stores the total number of replacements stored and to be executed. */
  uint16 m_numOfReplacements;

  /**Used as an escaping character. */
  byte m_escapeByte;

  /**Tells in what stage the algorithm is.
   * 0 - When initialized. Before hash table is resized.
   * 1 - During the analysation. After the hash table is resized.
   * 2 - Ready to decide replacements.
   */
  byte m_phase;

  const byte *m_data;
  
  /**Tells if the results of analysis and replacements are printed. */
  bool m_verbose;

  /**Tells if some rare symbols will be stolen. */
  const bool m_useEscaping;
};

} //namespace bwtc

#endif
