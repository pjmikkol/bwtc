/**
 * @file RunReplacer.hpp
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
 * Header for preprocessor which replaces the most often occurring runs of
 * same character.
 */
#ifndef RUN_REPLACER_HPP_
#define RUN_REPLACER_HPP_

#include "../globaldefs.hpp"
#include "FrequencyTable.hpp"
#include "../Utils.hpp"

#include <boost/static_assert.hpp>

#include <map>
#include <utility>
#include <vector>

namespace bwtc {

namespace RunReplacerConsts {
static const size_t s_maxLengthOfSequence = 1 << 15;
static const size_t s_logMaxLengthOfSequence = 15;
}

struct Runs {
  Runs(size_t freq, uint16 len, byte sym);
  Runs(const Runs& r);
  Runs& operator=(const Runs& r);
  
  size_t frequency;
  uint32 length;
  byte symbol;

  /**Does length have range big enough for storing the maximal runs? */
  BOOST_STATIC_ASSERT(RunReplacerConsts::s_maxLengthOfSequence <= (1 << 16) - 1);
};

bool operator<(const Runs& r1, const Runs& r2);

struct RunReplacement {
  RunReplacement(int next, size_t len, byte symbol);

  uint32 length;
  int16 nextElement;
  byte replacementSymbol;
  /**Does length have range big enough for storing the maximal runs? */
  BOOST_STATIC_ASSERT(RunReplacerConsts::s_maxLengthOfSequence <= (1 << 16) - 1);
};

/**This is used when constructing the RunReplacementTable. In that phase
 * nextElement stores actually the symbol of run represented by this
 * RunReplacement. */
bool operator<(const RunReplacement& r1, const RunReplacement& r2);

/**ReplacementTable stores the replacements for runs of same character.
 * It also stores the information about the characters to be escaped.
 * There exists four different options for character c:
 * * If runs of c aren't replaced and it isn't escaped, then
 *   m_listBegins[c] == -1 and m_escaped[c] == false
 * * If runs of c aren't replaced and it is escaped, then
 *   m_listBegins[c] == -1 and m_escaped[c] == true
 * * If runs of c are replaced and it isn't escaped, then
 *   m_listBegins[c] is the index of the replacement of longest run of c
 *   in array m_replacements. The second longest run is stored in index
 *   pointed by m_replacements[m_listBegins[c]].nextElement. Rest of the
 *   replacements are stored in same way. The value of last replacement's
 *   nextElement field is -1. In addition m_escaped[c] == false
 * * If runs of c are replaced and it is escaped, then
 *   same as above but m_escaped[c] == true.
 */
class RunReplacementTable {
 public:
  RunReplacementTable();

  /**Should be called when the replacements are stored in table.*/
  void prepare();

  void addReplacement(byte runSymbol, size_t length, byte replacementSymbol);
  
  void addEscaping(byte symbol);    
  
  bool noRunReplacements(byte symbol) const;

  bool isEscaped(byte symbol) const;

  int listPointer(byte symbol) const;

  const RunReplacement& replacement(int pointer) const;
  
  size_t size() const;

  bool escapingInUse() const;

  void printReplacements() const;

  void replacements(std::vector<std::pair<byte, RunReplacement> >& runReplacements) const;
  
 private:
  std::vector<RunReplacement> m_replacements;
  int16 m_listBegins[256];
  bool m_escaped[256];

};


class RunReplacer {
 public:
  RunReplacer();
  RunReplacer(bool verbose);
  RunReplacer(const RunReplacer& rr);
  ~RunReplacer();


  void analyseData(const byte *data, size_t length, bool reset=true);

  inline void analyseData(byte next);
  
  void beginAnalysing(byte first, bool reset);

  void finishAnalysation();

  size_t decideReplacements();

  size_t writeHeader(byte *to) const;

  size_t writeReplacedVersion(const byte *src, size_t length, byte *dst) const;

  void findReplaceableRuns(std::vector<Runs>& replaceableRuns,
                           FrequencyTable& freqs) const;

  size_t findEscapeIndex(FrequencyTable& freqs, size_t freeSymbols,
                         const std::vector<Runs>& replaceableRuns) const;
  
 private:
  void resetAnalyseData();
  inline void updateRunFrequency(std::pair<int, byte> run);
  size_t writeRunReplacement(byte runSymbol, int runLength, byte *dst) const;

  RunReplacer& operator=(const RunReplacer&);

  /**Stores the replacements after analysation phase. */
  RunReplacementTable m_replacements;

  /**Stores the frequencies of bytes. */
  size_t m_frequencies[256];

  /**Stores the frequencies of runs. The value stored in m_runFreqs[c][l],
   * is the amount of how many times the run of c of the length (1 << (l+1))
   * has been appeared. */
  std::vector<size_t> m_runFreqs[256];

  /**Stores the total number of replacements stored and to be executed. */
  uint16 m_numOfReplacements;

  /**Used in analysation phase. */
  std::pair<int, byte> m_prevRun;

  /**Used as an escaping character. */
  byte m_escapeByte;

  /**Tells if we have started the analysation. */
  bool m_analysationStarted;

  /**Tells if the results of analysis and replacements are printed. */
  bool m_verbose;
};

/**Class for storing the counts of runs, picking the most profitable run and
 * updating the frequencies of runs. Essentially this is heap with additional
 * keys, which provide fast search.
 */
class SequenceHeap {
 public:
  SequenceHeap();

  size_t maxUtility() const;

  bool empty() const;

  Runs removeMax();

  void prepare();

  void addRun(const Runs& run);
  
 private:
  void decreaseFrequency(int index, size_t value);

  void remove(int index);

  void initLocations();

  void heapify(int index);

  void buildMaxHeap();

  std::vector<Runs> m_runs;
  std::map<uint16, uint16> m_runLocations[256];
  int m_last;
  
  /** Is uint16 enough for the target of map? */
  
  BOOST_STATIC_ASSERT(RunReplacerConsts::s_logMaxLengthOfSequence*256 < (1 << 16) - 1);
  /** Is uint16 enough for the domain of map? */
  BOOST_STATIC_ASSERT(RunReplacerConsts::s_maxLengthOfSequence <= (1 << 16) - 1);
};

} //namespace bwtc

#endif
