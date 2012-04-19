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
#include "FrequencyTable.hpp"

#include <cassert>
#include <vector>
#include <utility>
#include <map>

namespace bwtc {

namespace long_sequences {

struct Sequence {
  Sequence(uint32 c, uint32 l, uint32 n, uint32 p)
      : count(c), length(l), name(n), samplePosition(p) {}

  bool operator<(const Sequence& s) const {
    int diff = count*(length-1) - s.count*(s.length-1);
    if(diff != 0) return diff < 0;
    else return name < s.name;
  }

  uint32 count;
  uint32 length;
  uint32 name;
  uint32 samplePosition;
};

} //namespace long_sequences

class SequenceReplacer {
 public:
  SequenceReplacer();
  SequenceReplacer(bool verbose);
  SequenceReplacer(const SequenceReplacer& sr);
  ~SequenceReplacer();

  void analyseData(const byte *data, size_t length, bool reset=true);

  void resetAnalyseData();

  void finishAnalysation();

  uint32 decideReplacements();

  size_t writeHeader(byte *to) const;

  size_t writeReplacedVersion(const byte *src, size_t length, byte *dst) const;

  struct bucket_struct {
    bucket_struct()
        : name(0), position(0xffffffff), positionInPos(0xffffffff) {}
    bucket_struct(uint32 n, uint32 p, uint32 pp)
        : name(n), position(p), positionInPos(pp) {}
    uint32 name;
    uint32 position;
    uint32 positionInPos;
  };

 private:
  void resizeAndInitTable(size_t preference);
  uint64 initHash(const byte* data) const;
  void initHashConstant();
  void scanAndStore();
  void calculateFrequencies(uint32 begin, uint32 end, uint32 *f) const;
  void sortIntoBuckets();
  void sortSubBucket(int begin, int end, bool positionsUnordered);
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
  void expandStringsInBucket(uint32 begin, uint32 end);
  std::pair<uint32, uint32> findLeftLimit(uint32 sequence, uint32 offset) const;
  std::pair<uint32, uint32> findRightLimit(uint32 sequence, uint32 offset) const;
  uint32 expandToLeft(const std::vector<uint32>& elements,
                      std::vector<std::pair<uint32, uint32> >& leftLimit);
  uint32 expandToRight(const std::vector<uint32>& elements,
                      std::vector<std::pair<uint32, uint32> >& rightLimit, uint32 leftExp);
  void removeOverlappingSequences(const std::vector<uint32>& elements, uint32 leftOffset,
                                  uint32 lengthOfSequence);
  void validateCorrectOrder(uint32 begin, uint32 end);
  void deleteRemovedAndTakeSamples(std::vector<long_sequences::Sequence>& replaceables);
  uint32 findReplaceableSequences(const std::vector<long_sequences::Sequence>& replaceables,
                                  FrequencyTable& freqTable, uint32 maxReps) const;
  uint32 findEscapeIndex(FrequencyTable& freqTable, uint32 freeSymbols,
                         std::vector<long_sequences::Sequence>& replaceables,
                         uint32 candidates) const;
  uint32 writeAndPackInteger(byte *to, uint32 length) const;
  uint32 writeWithEscaping(const byte* begin, const byte* end, byte* dst) const;

  
  static const uint64 s_hashConstant = 37;
  static const uint32 s_errorVal = 0xffffffff;
  static const uint32 s_defaultWindowSize = 32;
  static const uint32 s_maxExpandableBucket = 10;
  //static const int s_insertionSortLimit = 10;

  /**Stores the frequencies of bytes. */
  uint32 m_frequencies[256];

  bool m_isEscaped[256];

  /**Stores the hash values, counters and lengths of replaceable sequences.
   * When scanning, the first one is counter and second is h* */
  std::vector<std::pair<uint32, uint32> > m_hashValues;

  std::vector<std::pair<uint32, uint32> > m_sequences;

  std::vector<bucket_struct> m_buckets;

  std::map<uint32, const byte*> m_samples;

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

};

} //namespace bwtc

#endif
