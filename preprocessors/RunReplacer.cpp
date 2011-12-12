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
 * Implementation of RunReplacer.
 */
#include <algorithm>
#include <utility>

#include "RunReplacer.hpp"
#include "../Utils.hpp"

namespace bwtc {

Runs::Runs(size_t freq, uint16 len, byte sym)
  : frequency(freq), length(len), symbol(sym) {}

Runs::Runs(const Runs& r)
    : frequency(r.frequency), length(r.length), symbol(r.symbol) {}

Runs& Runs::operator=(const Runs& r) {
  frequency = r.frequency;
  length = r.length;
  symbol = r.symbol;
  return *this;
}

bool operator<(const Runs& r1, const Runs& r2) {
  return (r1.length - 1)*r1.frequency < (r2.length - 1)*r2.frequency;
}

RunReplacement::RunReplacement(int next, size_t len, byte symbol)
    : length(len), nextElement(next), replacementSymbol(symbol) {}

bool operator<(const RunReplacement& r1, const RunReplacement& r2) {
  return r1.nextElement < r2.nextElement ||
      (r1.nextElement == r2.nextElement && r1.length < r2.length);
}


SequenceHeap::SequenceHeap() : m_last(-1) {}

void SequenceHeap::addRun(const Runs& run) {
  m_runs.push_back(run);
  ++m_last;
  assert(m_last == ((int)m_runs.size()) - 1);
}

size_t SequenceHeap::maxUtility() const {
  const Runs& max = m_runs[0];
  return (max.length - 1)*max.frequency;
}

bool SequenceHeap::empty() const { return m_last < 0; }

void SequenceHeap::prepare() {
  initLocations();
  assert(m_last == ((int)m_runs.size())-1);
  buildMaxHeap();
}

Runs SequenceHeap::removeMax() {
  Runs max = m_runs[0];
  for(std::map<uint16, uint16>::iterator it = m_runLocations[max.symbol].begin();
      it != m_runLocations[max.symbol].end(); ++it)
  {
    if(it->first >= max.length) remove(it->second);
  }
  return max;
}

void SequenceHeap::remove(int index) {
  if(index > m_last) return;
  m_runLocations[m_runs[index].symbol][m_runs[index].length] = m_last;
  m_runLocations[m_runs[m_last].symbol][m_runs[m_last].length] = index;
  std::swap(m_runs[index], m_runs[m_last]);
  --m_last;
  heapify(index);
}

void SequenceHeap::decreaseFrequency(int index, size_t value) {
  if(index > m_last) return;
  //assert(m_runs[index].frequency >= value);
  //m_runs[index].frequency -= value;
  if(m_runs[index].frequency >= value) m_runs[index].frequency -= value;
  else m_runs[index].frequency = 0;
  heapify(index);
}

void SequenceHeap::initLocations() {
  for(size_t i = 0; i < m_runs.size(); ++i) {
    const Runs& r = m_runs[i];
    m_runLocations[r.symbol].insert(std::make_pair(r.length, i));
    assert(r.length > 0);
  }
}
#define parent(x) (((x)-1) >> 1)
#define left(x) (2*(x) + 1)
#define right(x) (2*(x) + 2)


void SequenceHeap::heapify(int i) {
  int l = left(i), r = right(i);
  while(r <= m_last) {
    int largest = (m_runs[l] < m_runs[r])? r : l;
    if(m_runs[i] < m_runs[largest]) {
      const Runs& lr = m_runs[largest];
      const Runs& ir = m_runs[i];
      assert(m_runLocations[lr.symbol].count(lr.length) > 0);
      assert(m_runLocations[ir.symbol].count(ir.length) > 0);
      std::swap(m_runLocations[lr.symbol][lr.length],
                m_runLocations[ir.symbol][ir.length]);
      std::swap(m_runs[largest], m_runs[i]);
      i = largest;
      l = left(i), r = right(i);
    } else {
      return;
    }
  }
  if(l == m_last && m_runs[i] < m_runs[l]) {
      std::swap(m_runLocations[m_runs[l].symbol][m_runs[l].length],
                m_runLocations[m_runs[i].symbol][m_runs[i].length]);
      std::swap(m_runs[i], m_runs[l]);
  }
}

void SequenceHeap::buildMaxHeap() {
  for(int i = parent(m_last); i >= 0; --i) heapify(i);
}

#undef parent
#undef left
#undef right

// End of SequenceHeap's implementation


// Implementation of RunReplacementTable

RunReplacementTable::RunReplacementTable() {
  std::fill(m_escaped, m_escaped + 256, false);
  std::fill(m_listBegins, m_listBegins + 256, -1);
}

void RunReplacementTable::prepare() {
  std::sort(m_replacements.begin(), m_replacements.end());
  for(size_t i = 0; i < m_replacements.size(); ++i) {
    byte symbol = (byte)m_replacements[i].nextElement;
    int oldPointer = m_listBegins[symbol];
    m_listBegins[symbol] = i;
    m_replacements[i].nextElement = oldPointer;
  }
}

void RunReplacementTable::addReplacement(byte runSymbol, size_t length,
                                         byte replacementSymbol)
{
  m_replacements.push_back(RunReplacement(runSymbol,length,replacementSymbol));
}

size_t RunReplacementTable::size() const { return m_replacements.size(); }

bool RunReplacementTable::noRunReplacements(byte symbol) const {
  return m_listBegins[symbol] == -1;
}

bool RunReplacementTable::isEscaped(byte symbol) const {
  return m_escaped[symbol];
}

void RunReplacementTable::addEscaping(byte symbol) {
  m_escaped[symbol] = true;
}

int RunReplacementTable::listPointer(byte symbol) const {
  return m_listBegins[symbol];
}

const RunReplacement& RunReplacementTable::replacement(int pointer) const {
  assert(pointer >= 0);
  assert(pointer < (int)m_replacements.size());
  return m_replacements[pointer];
}

bool RunReplacementTable::escapingInUse() const {
  for(size_t i = 0; i < 256; ++i) {
    if(m_escaped[i]) return true;
  }
  return false;
}

void RunReplacementTable::printReplacements() const {
  for(size_t s = 0; s < 255; ++s) {
    for(int i = m_listBegins[s]; i != -1; i = m_replacements[i].nextElement) 
      std::cout << (int) m_replacements[i].replacementSymbol << " -> "
                << m_replacements[i].length << "x" << s << std::endl;
  }
}

void RunReplacementTable::
replacements(std::vector<std::pair<byte, RunReplacement> >& runReplacements) const {
  for(size_t s = 0; s < 255; ++s) {
    for(int i = m_listBegins[s]; i != -1; i = m_replacements[i].nextElement) 
      runReplacements.push_back(std::make_pair((byte) s, m_replacements[i]));
  }
}

// End of RunReplacementTable's implementation

RunReplacer::RunReplacer()
    : m_numOfReplacements(0), m_prevRun(std::make_pair(0, 0)),
      m_analysationStarted(false), m_verbose(false)
{
  size_t size = RunReplacerConsts::s_logMaxLengthOfSequence - 1;
  for(size_t i = 0; i < 256; ++i) 
    m_runFreqs[i].resize(size);
  resetAnalyseData();
}

RunReplacer::RunReplacer(bool verbose)
    : m_numOfReplacements(0), m_prevRun(std::make_pair(0, 0)),
      m_analysationStarted(false), m_verbose(verbose)
{
  size_t size = RunReplacerConsts::s_logMaxLengthOfSequence - 1;
  for(size_t i = 0; i < 256; ++i) 
    m_runFreqs[i].resize(size);
  resetAnalyseData();
}

RunReplacer::RunReplacer(const RunReplacer& rr)
    : m_numOfReplacements(rr.m_numOfReplacements), m_prevRun(rr.m_prevRun),
      m_analysationStarted(rr.m_analysationStarted), m_verbose(rr.m_verbose)
{
  std::copy(rr.m_frequencies, rr.m_frequencies + 256, m_frequencies);
  for(size_t i = 0; i < 256; ++i) {
    m_runFreqs[i] = rr.m_runFreqs[i];
  }
}

RunReplacer::~RunReplacer() {}

void RunReplacer::resetAnalyseData() {
  std::fill(m_frequencies, m_frequencies + 256, 0);
  for(size_t i = 0; i < 256; ++i) 
    std::fill(m_runFreqs[i].begin(), m_runFreqs[i].end(), 0);
  m_prevRun = std::make_pair(0, 0);
}

void RunReplacer::beginAnalysing(byte first, bool reset) {
  if(reset) {
    resetAnalyseData();
  }
  m_analysationStarted = true;
  m_prevRun = std::make_pair(1, first);
  ++m_frequencies[first];
}

void RunReplacer::updateRunFrequency(std::pair<uint16, byte> run) {
  assert(run.first <= RunReplacerConsts::s_maxLengthOfSequence);
  assert(run.first > 1);
  run.first &= 0xfffffffe;
  while(run.first) {
    int logLongest = utils::logFloor((unsigned)run.first);
    ++m_runFreqs[run.second][logLongest - 1];
    run.first ^= (1 << logLongest);
  }
}

void RunReplacer::findReplaceableRuns(std::vector<Runs>& replaceableRuns,
                                      FrequencyTable& freqs) const
{
  assert(replaceableRuns.size() == 0);
  SequenceHeap seqHeap;
  for(size_t i = 0; i < 256; ++i) {
    for(size_t j = 0; j < m_runFreqs[i].size(); ++j) {
      if(m_runFreqs[i][j] > 0) {
        seqHeap.addRun(Runs(m_runFreqs[i][j], 1 << (j+1), (byte)i));
      }
    }
  }
  seqHeap.prepare();
  size_t currentSymbol = 0;
  while(!seqHeap.empty()) {
    Runs best = seqHeap.removeMax();
    size_t score = best.length*best.frequency;
    freqs.decrease(best.symbol, score);
    if(freqs.getFrequency(currentSymbol) + 3 >= score - best.frequency) {
      freqs.increase(best.symbol, score);
      break;
    }
    replaceableRuns.push_back(best);
    ++currentSymbol;
  }
}

size_t RunReplacer::findEscapeIndex(FrequencyTable& freqs, size_t freeSymbols,
                                    const std::vector<Runs>& replaceableRuns) const
{
  if(replaceableRuns.size() <= freeSymbols) return freeSymbols;
  assert(replaceableRuns.size() <= 254);
  int64 utility = 0;
  size_t i = freeSymbols;
  for(; i < replaceableRuns.size(); ++i) {
    const Runs& r = replaceableRuns[i];
    utility += ((r.length - 1)*r.frequency - freqs.getFrequency(i) - 3);
  }
  while(utility <= ((int64)freqs.getFrequency(i)) && i > freeSymbols) {
    --i;
    const Runs& r = replaceableRuns[i];
    int64 score = (r.length - 1)*r.frequency;
    freqs.increase(r.symbol, score + r.frequency);
    //freqs.increase(r.symbol, score);
    utility -= (score - freqs.getFrequency(i) - 3);
  }
  return i;
}

size_t RunReplacer::writeHeader(byte *to) const {
  assert(m_replacements.size() == m_numOfReplacements);
  if(m_numOfReplacements == 0) {
    to[0] = to[1] = 0;
    return 2;
  }

  std::vector<std::pair<byte, RunReplacement> > replacements;
  m_replacements.replacements(replacements);
  assert(replacements.size() == m_numOfReplacements);
  
  size_t pos = 0;
  size_t pairs = m_numOfReplacements & 0xfffffffe;
  byte prev = m_escapeByte;
  for(size_t i = 0; i < pairs; i += 2) {
    to[pos++] = replacements[i].second.replacementSymbol;
    byte lengths = (utils::logFloor((uint32)replacements[i].second.length) << 4) |
        (utils::logFloor((uint32)replacements[i+1].second.length));
    to[pos++] = lengths;
    to[pos++] = replacements[i].first;
    to[pos++] = replacements[i+1].second.replacementSymbol;
    to[pos++] = replacements[i+1].first;
  }

  if(m_numOfReplacements != pairs) {
    prev = replacements.back().second.replacementSymbol;
    to[pos++] = prev;
    to[pos++] = utils::logFloor((uint32)replacements.back().second.length) << 4;
    to[pos++] = replacements.back().first;
    to[pos++] = (m_replacements.escapingInUse())?m_escapeByte:prev;
  } else {
    to[pos++] = (m_replacements.escapingInUse())?m_escapeByte:prev;
    to[pos++] = 0;
  }
  return pos;
}

size_t RunReplacer::
writeRunReplacement(byte runSymbol, int runLength, byte *dst) const {
  size_t j = 0;
  int pointer = m_replacements.listPointer(runSymbol);
  while(pointer != -1 && runLength > 0) {
    const RunReplacement& rr = m_replacements.replacement(pointer);
    if((int)rr.length <= runLength) {
      // TODO: optimize (rr.length == 2^k, for some k)
      size_t times = runLength/rr.length;
      std::fill(dst +j, dst + j + times, rr.replacementSymbol);
      runLength -= times*rr.length;
      j += times;
    }
    pointer = rr.nextElement;
  }
  if(runLength > 0) {
    if(m_replacements.isEscaped(runSymbol)) {
      for(int i = 0; i < runLength; ++i) {
        dst[j++] = m_escapeByte; dst[j++] = runSymbol;
      }
    } else {
      std::fill(dst + j, dst + j + runLength, runSymbol);
      j += runLength;
    }
  }
  return j;
}


size_t RunReplacer::
writeReplacedVersion(const byte *src, size_t length, byte *dst) const {
  size_t j = 0;
  byte prev = src[0];
  for(size_t i = 1; i < length; ++i) {
    if(prev != src[i]) {
      if(m_replacements.isEscaped(prev)) dst[j++] = m_escapeByte;
      dst[j++] = prev;
    } else {
      int runLength = 2;
      ++i;
      while(i < length && prev == src[i]) {++i; ++runLength; }
      j += writeRunReplacement(prev, runLength, dst + j);
    }
    prev = src[i];
  }
  if(prev != src[length - 2]) {
    if(m_replacements.isEscaped(prev)) dst[j++] = m_escapeByte;
    dst[j++] = prev;
  }
  return j;
}

size_t RunReplacer::decideReplacements() {
  FrequencyTable freqTable(m_frequencies);
  size_t freeSymbols = 0;
  while(freqTable.getFrequency(freeSymbols) == 0) ++freeSymbols;

  std::vector<Runs> replaceableRuns;
  findReplaceableRuns(replaceableRuns, freqTable);

  size_t escapeIndex = (replaceableRuns.size() <= freeSymbols)?freeSymbols:
      findEscapeIndex(freqTable, freeSymbols, replaceableRuns);

  m_numOfReplacements = (escapeIndex > freeSymbols)?
      escapeIndex:std::min(freeSymbols, replaceableRuns.size());

  if(m_verbose) {
    std::clog << "Replacing " << m_numOfReplacements << " runs. ";
    if(m_numOfReplacements > freeSymbols)
      std::clog << "Made " << (escapeIndex - freeSymbols + 1)
                << " symbols free." << std::endl;
    else
      std::clog << "No symbols made free." << std::endl;
  }
  if(escapeIndex > freeSymbols) {
    for(size_t i = freeSymbols; i <= escapeIndex; ++i)
      m_replacements.addEscaping(freqTable.getKey(i));
  }
  for(size_t i = 0; i < m_numOfReplacements; ++i) {
    const Runs& r = replaceableRuns[i];
    m_replacements.addReplacement(r.symbol, r.length, freqTable.getKey(i));
  }
  m_replacements.prepare();
  return m_replacements.size();
}

void RunReplacer::analyseData(const byte* data, size_t length, bool reset) {
  size_t i = 0;
  if(!m_analysationStarted) {
    beginAnalysing(data[i++], reset);
  }
  for(;i < length; ++i) {
    analyseData(data[i]);
  }
}

void RunReplacer::analyseData(byte next) {
  assert(m_analysationStarted);
  if (next == m_prevRun.second &&
      m_prevRun.first < (int)RunReplacerConsts::s_maxLengthOfSequence)
  {
    ++m_prevRun.first;
  } else {
    if(m_prevRun.first > 1) {
      updateRunFrequency(m_prevRun);
    }
    m_prevRun = std::make_pair(1, next);
  }
  ++m_frequencies[next];
}

void RunReplacer::finishAnalysation() {
  if(m_prevRun.first > 1) {
      updateRunFrequency(m_prevRun);
  }
  m_prevRun = std::make_pair(0,0);
}

} //namespace bwtc
