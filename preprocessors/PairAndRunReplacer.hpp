/**
 * @file PairAndRunReplacer.hpp
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
 * Header for preprocessor which replaces the most often occurring pairs
 * and runs of the same character.
 */
#ifndef PAIR_AND_RUN_REPLACER_HPP_
#define PAIR_AND_RUN_REPLACER_HPP_

#include "../globaldefs.hpp"
#include "FrequencyTable.hpp"
#include "RunReplacer.hpp"

#include <vector>
#include <utility>

namespace bwtc {

namespace pairs_and_runs {
struct Replacement {
  Replacement(size_t freq, uint32 len, uint16 repl, bool isPair);
  Replacement(const Replacement& r);
  Replacement& operator=(const Replacement& r);
  
  size_t frequency;
  /**Length of run, if the replacement is for run. */
  uint32 length;
  /**Pair or the symbol to be replaced. */
  uint16 replaceable;
  bool pair;
};

class PairAndRunReplacer {
 public:
  PairAndRunReplacer(bool useEscaping);
  PairAndRunReplacer(bool useEscaping, bool verbose);
  PairAndRunReplacer(const PairAndRunReplacer& pr);
  ~PairAndRunReplacer();

  void analyseData(const byte *data, size_t length, bool reset=true);

  inline void analyseData(byte next);
  
  void beginAnalysing(byte first, bool reset);

  void finishAnalysation();

  void findReplaceablePairsAndRuns(
      std::vector<std::pair<size_t, uint16> >& pairs,
      std::vector<Replacement>& replaceables, FrequencyTable& freqs,
      size_t maxReplaceables) const;

  size_t findEscapeIndex(FrequencyTable& freqs, size_t freeSymbols,
                         std::vector<Replacement>& replaceables);
  
  size_t decideReplacements();
  
  size_t writeHeader(byte *dst) const;

  size_t writeReplacedVersion(const byte *src, size_t length, byte *dst) const;

  size_t writeRunReplacement(byte runSymbol, int runLength, byte *dst) const;

  inline void updateRunFrequency(std::pair<int, byte> run);
  
 private:
  void resetAnalyseData();
  PairAndRunReplacer& operator=(const PairAndRunReplacer&);

  /**Stores the frequencies of bytes. */
  uint32 m_frequencies[256];

  /**Stores the frequencies of runs. The value stored in m_runFreqs[c][l],
   * is the amount of how many times the run of c of the length (1 << (l+1))
   * has been appeared. */
  std::vector<size_t> m_runFreqs[256];

  /**Used for storing the counters when analyzing data. */
  size_t m_pairFrequencies[1 << 16];

  /**Stores the pair replacements after analysation phase. */
  RunReplacementTable m_runReplacements;


  /** In replacement writing phase the replacements are stored in this. */
  byte m_pairReplacements[1 << 16];

  /**Used in the analysation phase. */
  uint16 m_prev;

  /**Used in analysation phase. */
  std::pair<int, byte> m_prevRun;

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

  const bool m_useEscaping;
};

} //namespace pairs_and_runs
} //namespace bwtc

#endif
