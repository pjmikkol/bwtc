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

#include "../globaldefs.hpp"
#include "FrequencyTable.hpp"

#include <cassert>
#include <vector>
#include <utility>

namespace bwtc {

class PairReplacer {
 public:
  PairReplacer();
  PairReplacer(bool verbose);
  PairReplacer(const PairReplacer& pr);
  ~PairReplacer();

  void analyseData(const byte *data, size_t length, bool reset=true);

  inline void analyseData(byte next) {
      assert(m_analysationStarted);
      m_prev = (m_prev << 8) | next;
      
      ++m_pairFrequencies[m_prev];
      ++m_frequencies[next];
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
                            FrequencyTable& freqs) const;


  size_t findEscapeIndex(FrequencyTable& freqs, size_t freeSymbols,
                     std::vector<std::pair<size_t, uint16> >& suitablePairs);

  size_t writeHeader(byte *to) const;

  size_t writeReplacedVersion(const byte *src, size_t length, byte *dst) const;

  template <typename PreProc>
  size_t writeAndAnalyseHeader(byte *to, PreProc& pp) const {
    pp.beginAnalysing(true);
    if(m_numOfReplacements == 0)  {
      to[0] = to[1] = to[2] = 0;
      pp.analyseData(0); pp.analyseData(0); pp.analyseData(0);
      return 3;
    }
    size_t pos = 0;
    byte prevValue = m_escapeByte;
    for(size_t i = 0; i < (1 << 16); ++i) {
      byte val = m_replacements[i];
      if(val != m_commonByte && val != m_escapeByte) {
        prevValue = val;
        to[pos] = prevValue;
        pp.analyseData(to[pos++]);
        to[pos] = (i >> 8) & 0xff;
        pp.analyseData(to[pos++]);
        to[pos] = i & 0xff;
        pp.analyseData(to[pos++]);
      }
    }
    to[pos] = prevValue;
    pp.analyseData(to[pos++]);
    to[pos] = (m_escapeByte != m_commonByte)?m_escapeByte:prevValue;
    pp.analyseData(to[pos++]);
    return pos;    
  }

  template <typename PreProc>
  size_t writeAndAnalyseReplacedVersion(
      const byte *src, size_t length, byte *dst, PreProc& pp) const
  {
    size_t j = 0; /* Is used for indexing the target. */
    uint16 pair = src[0];
    size_t i = 1;
    while(true) {
      pair = (pair << 8) | src[i];
      byte replValue = m_replacements[pair];
      if(replValue == m_commonByte) {
        dst[j] = src[i-1];
        pp.analyseData(dst[j++]);
      } else if(replValue != m_escapeByte) {
        dst[j] = replValue;
        pp.analyseData(dst[j++]);
        if(i == length - 1) break;
        pair = src[++i];
      } else {
        dst[j] = m_escapeByte;
        pp.analyseData(dst[j++]);
        dst[j] = src[i-1];
        pp.analyseData(dst[j++]);
      }
      if(i == length - 1) {
        pair = (src[i] << 8) | 0;
        if(m_replacements[pair] == m_escapeByte && m_escapeByte != m_commonByte) {
          dst[j] = m_escapeByte;
          pp.analyseData(dst[j++]);
        }
        dst[j] = src[i];
        pp.analyseData(dst[j++]);
        break;
      }
      ++i;
    }
    return j;
  }

 private:
  PairReplacer& operator=(const PairReplacer&);
  void constructReplacementTable(
      const std::vector<std::pair<size_t, uint16> >& pairs,
      const FrequencyTable& freqTable, size_t freeSymbols);

  /**Stores the frequencies of bytes. */
  size_t m_frequencies[256];

  /**Used for storing the counters when analyzing data. */
  size_t m_pairFrequencies[1 << 16];

  /** In replacement writing phase the replacements are stored in this. */
  byte m_replacements[1 << 16];

  /**Used in the analysation phase. */
  uint16 m_prev;

  /** Stores the total number of replacements stored and to be executed. */
  uint16 m_numOfReplacements;

  /**Used as an escaping character. */
  byte m_escapeByte;

  /**Used in writing phase as an indicator for "no replacement". */
  byte m_commonByte;
  
  /**Tells if the analysation is started. Is used for checks in
   * analyseData(byte* data, size_t length)-function: to be able to
   * know if there is something significant stored to m_prev.
   * analyseData(byte next) asserts that analysation has been started.*/
  bool m_analysationStarted;

  /**Tells if the results of analysis and replacements are printed. */
  bool m_verbose;
};

} //namespace bwtc

#endif
