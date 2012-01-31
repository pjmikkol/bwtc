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
#include <map>

namespace bwtc {

namespace long_sequences {

template <typename K, typename V>
class MaxHeap {
 public:
  MaxHeap() {}

  void insert(K key, V value) {
    if(m_positions.find(key) != m_positions.end()) {
      // check the heap-condition:
      m_positions[key] = value;
    } else {
      m_heap.push_back(std::make_pair(key, value));
    }
  }

 private:
  std::vector<std::pair<K, V> > m_heap;
  std::map<K, size_t> m_positions;
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
  void scanAndStore(const byte* data, size_t length);
  void calculateFrequencies(const byte* data, uint32 begin, uint32 end);
  void sortIntoBuckets();
  void sortSubBucket(int begin, int end, const byte* data);
  int strCmp(uint32 pos1, uint32 pos2, const byte* data) const;
  void sortPositions(int begin, int end);
  void insertionSort(int begin, int end, const byte* data);
  void sortAndMarkBuckets(const byte* data);
  /**Names hash values and removes useless strings. Returns number of
   * separate string.*/
  uint32 nameHashValues();
  void nameRange(uint32 begin, uint32 end, uint32 name);
  bool validatePhase2(const byte* data) const;
  bool validateRange(uint32 begin, uint32 end, const byte* data) const;
  void removeDeletedSequences();

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

  uint32 m_windowSize;
  
  /** Stores the total number of replacements stored and to be executed. */
  uint16 m_numOfReplacements;

  /**Used as an escaping character. */
  byte m_escapeByte;

  /**Tells in what stage the algorithm is.
   * 0 - When initialized. Before hash table is resized.
   * 1 - During the analysation. After the hash table is resized.
   * 2 - Bucket sorting phase.
   * 3 - Bucket analysation phase.
   * 4 - Deletion of the overlapping strings.
   * 5 - Final analysation of replacements.
   * 6 - Writing of replaced string.
   */
  byte m_phase;

  /**Tells if the results of analysis and replacements are printed. */
  bool m_verbose;

  /**Tells if some rare symbols will be stolen. */
  const bool m_useEscaping;
};

} //namespace bwtc

#endif
