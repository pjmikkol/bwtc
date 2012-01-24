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

  static const uint64 s_hashConstant = 37;
  
  /**Stores the frequencies of bytes. */
  size_t m_frequencies[256];

  /**Stores the hash values, counters and lengths of replaceable sequences. */
  std::vector<std::pair<uint32, uint32> > m_hashValues;

  uint64 m_runningHash;

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
