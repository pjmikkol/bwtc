/**
 * @file PairReplacer.hpp
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
 * Header for preprocessor which replaces the most often occurring pairs.
 */
#ifndef PAIR_REPLACER_HPP_
#define PAIR_REPLACER_HPP_

#include <cassert>
#include <vector>
#include <utility>

#include "../globaldefs.hpp"
#include "FrequencyTable.hpp"
#include "Grammar.hpp"

namespace bwtc {

class PairReplacer {
 public:
  PairReplacer(Grammar& grammar);
  PairReplacer(Grammar& grammar, bool verbose);
  ~PairReplacer();

  void analyseData(const byte *data, size_t length, bool reset=true);

  inline void analyseData(byte next) {
    assert(m_analysationStarted);
    assert(m_pairFrequencies);
    assert(m_frequencies);
    m_prev = (m_prev << 8) | next;
    
    ++m_pairFrequencies[m_prev];
    ++m_frequencies[next];
  }

  inline void analyseData0(byte next) {
    ++m_frequencies[next];
    uint16 prev = m_prev;
    m_prev = (m_prev << 8) | next;

    int pXc = prev ^ m_prev;
    // compute min(pXc,1)
    int min = 1 + ((pXc - 1) & ((pXc - 1) >> (sizeof(int) * 8 - 1)));
    
    m_pairFrequencies[m_prev] += min;
  }

  void beginAnalysing(bool reset);

  void beginAnalysing(byte first, bool reset);

  void resetAnalyseData();

  void finishAnalysation();

  size_t decideReplacements();
  
  static void makePairList(std::vector<std::pair<size_t, uint16> >& pairs,
                           const size_t *pairFrequencies);

  void findReplaceablePairs(std::vector<std::pair<size_t, uint16> >& pairs,
                            std::vector<std::pair<size_t, uint16> >& replaceablePairs,
                            FrequencyTable& freqs, size_t maxReplacements,
                            uint32& variables, uint32& specials,
                            uint32& forFree) const;

  int64 findReplaceables(size_t start, 
                         const std::vector<std::pair<size_t, uint16> >& pairs,
                         std::vector<std::pair<size_t, uint16> >& replPairs,
                         FrequencyTable& freqs, size_t maxRepl,
                         uint32& variables, uint32& specials,
                         uint32& forFree) const;

  size_t writeReplacedVersion(const byte *src, size_t length, byte *dst) const;


 private:
  PairReplacer& operator=(const PairReplacer&);
  PairReplacer(const PairReplacer& pr);
  void constructReplacementTable(
      const std::vector<std::pair<size_t, uint16> >& pairs,
      const std::vector<byte>& freedSymbols,
      const std::vector<byte>& newSpecials,
      const std::vector<byte>& replacements);

  Grammar& m_grammar;

  /**Stores the frequencies of bytes. */
  size_t m_frequencies[256];

  /**Used for storing the counters when analyzing data. */
  size_t m_pairFrequencies[1 << 16];

  /** In replacement writing phase the replacements are stored in this.
   *  - if both of the bytes are m_commonByte nothing is to be done
   *  - if both are different than m_commonByte first character is to
   *    replaced with pair in that cell
   *  - if second byte is not common byte pair is replaced with second byte
   */
  uint16 m_replacements[1 << 16];

  /**Used in the analysation phase. */
  uint16 m_prev;

  /** Stores the total number of replacements stored and to be executed. */
  uint16 m_numOfReplacements;

  uint16 m_numOfFreedSymbols;
  uint16 m_numOfNewSpecials;
  
  /**Used in writing phase as an indicator for "no replacement". */
  byte m_commonByte;
  
  /**Tells if the analysation is started. Is used for checks in
   * analyseData(byte* data, size_t length)-function: to be able to
   * know if there is something significant stored to m_prev.
   * analyseData(byte next) asserts that analysation has been started.*/
  bool m_analysationStarted;

  /**Tells if the results of analysis and replacements are printed. */
  bool m_verbose;

  /**Number of different greedy searchs to be made when deciding
   * replaceable pairs.*/
  static const size_t s_greedyStarts = 5;
};

} //namespace bwtc

#endif
