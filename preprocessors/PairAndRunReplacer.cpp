/**
 * @file PairAndRunReplacer.cpp
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
 * Implementation of preprocessor which replaces the most often occurring pairs
 * and runs of the same character.
 */


#include "../globaldefs.hpp"
#include "PairAndRunReplacer.hpp"
#include "RunReplacer.hpp"
#include "PairReplacer.hpp"
#include "FrequencyTable.hpp"

#include <algorithm>
#include <vector>
#include <utility>

namespace bwtc {
namespace pairs_and_runs {

Replacement::Replacement(size_t freq, uint32 len, uint16 repl, bool isPair)
    : frequency(freq), length(len), replaceable(repl), pair(isPair) {}

Replacement::Replacement(const Replacement& r)
    : frequency(r.frequency), length(r.length), replaceable(r.replaceable),
      pair(r.pair) {}

Replacement& Replacement::operator=(const Replacement& r) {
  frequency = r.frequency;
  length = r.length;
  replaceable = r.replaceable;
  pair = r.pair;
  return *this;
}

PairAndRunReplacer::PairAndRunReplacer(bool useEscaping)
    : m_prev(0), m_prevRun(std::make_pair(0,0)), m_numOfReplacements(0),
      m_escapeByte(0), m_commonByte(0), m_analysationStarted(false),
      m_verbose(false), m_useEscaping(useEscaping) {
  size_t size = RunReplacerConsts::s_logMaxLengthOfSequence - 1;
  for(size_t i = 0; i < 256; ++i) 
    m_runFreqs[i].resize(size);
}

PairAndRunReplacer::PairAndRunReplacer(bool useEscaping, bool verbose)
    : m_prev(0), m_prevRun(std::make_pair(0,0)), m_numOfReplacements(0),
      m_escapeByte(0), m_commonByte(0), m_analysationStarted(false),
      m_verbose(verbose), m_useEscaping(useEscaping) {
  size_t size = RunReplacerConsts::s_logMaxLengthOfSequence - 1;
  for(size_t i = 0; i < 256; ++i) 
    m_runFreqs[i].resize(size);
}

PairAndRunReplacer::PairAndRunReplacer(const PairAndRunReplacer& pr)
    : m_runFreqs(pr.m_runFreqs), m_runReplacements(pr.m_runReplacements),
      m_prev(pr.m_prev), m_prevRun(pr.m_prevRun),
      m_numOfReplacements(pr.m_numOfReplacements), m_escapeByte(pr.m_escapeByte),
      m_commonByte(pr.m_commonByte), m_analysationStarted(pr.m_analysationStarted),
      m_verbose(pr.m_verbose), m_useEscaping(pr.m_useEscaping) {
  std::copy(pr.m_frequencies, pr.m_frequencies + 256, m_frequencies);
  std::copy(pr.m_pairFrequencies, pr.m_pairFrequencies + (1 << 16), m_pairFrequencies);
  for(size_t i = 0; i < 256; ++i) {
    m_runFreqs[i] = pr.m_runFreqs[i];
  }
}

PairAndRunReplacer::~PairAndRunReplacer() {}

void PairAndRunReplacer::resetAnalyseData() {
  std::fill(m_frequencies, m_frequencies + 256, 0);
  std::fill(m_pairFrequencies, m_pairFrequencies + (1 << 16), 0);
}

void PairAndRunReplacer::analyseData(const byte *data, size_t length, bool reset) {
  assert(length > 0);
  size_t i = 0;
  if(!m_analysationStarted) beginAnalysing(data[i++], reset);
  for(; i < length; ++i) {
    analyseData(data[i]);
  }
}

inline void PairAndRunReplacer::updateRunFrequency(std::pair<int, byte> run) {
  assert(run.first <= (int)RunReplacerConsts::s_maxLengthOfSequence);
  assert(run.first > 1);
  run.first &= 0xfffffffe;
  while(run.first) {
    int logLongest = utils::logFloor((unsigned)run.first);
    ++m_runFreqs[run.second][logLongest - 1];
    run.first ^= (1 << logLongest);
  }
}


inline void PairAndRunReplacer::analyseData(byte next) {
  assert(m_analysationStarted);
  m_prev = (m_prev << 8) | next;

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

  ++m_pairFrequencies[m_prev];
  ++m_frequencies[next];
}

void PairAndRunReplacer::finishAnalysation() {
  if(m_prevRun.first > 1) {
      updateRunFrequency(m_prevRun);
  }
  m_prevRun = std::make_pair(0,0);
}

void PairAndRunReplacer::beginAnalysing(byte first, bool reset) {
  if(reset) resetAnalyseData();
  m_analysationStarted = true;
  m_prevRun = std::make_pair(1, first);
  m_prev = first;
  ++m_frequencies[first];
}

void PairAndRunReplacer::findReplaceablePairsAndRuns(
    std::vector<std::pair<size_t, uint16> >& pairs,
    std::vector<Replacement>& replaceables, FrequencyTable& freqs,
    size_t maxReplaceables) const {
  bool usedFst[256] = {false}, usedSnd[256] = {false}, usedRun[256] = {false};
  SequenceHeap seqHeap;
  for(size_t i = 0; i < 256; ++i) {
    for(size_t j = 0; j < m_runFreqs[i].size(); ++j) {
      if(m_runFreqs[i][j] > 0) {
        seqHeap.addRun(Runs(m_runFreqs[i][j], 1 << (j+1), (byte)i));
      }
    }
  }
  seqHeap.prepare();
  std::sort(pairs.rbegin(), pairs.rend());

  size_t currentSymbol = 0, currentPair = 0;

  while(currentSymbol < maxReplaceables) {
    bool isPair = true;
    size_t utility = 0;
    bool runsLeft = !seqHeap.empty(), pairsLeft = currentPair < pairs.size();
    if (runsLeft && pairsLeft) {
      size_t runUtility = seqHeap.maxUtility();
      size_t pairUtility = pairs[currentPair].first;
      if(runUtility < pairUtility) {
        isPair = true; utility = pairUtility;
      } else {
        isPair = false; utility = runUtility;
      }
    } else if(pairsLeft) {
      isPair = true;
      utility = pairs[currentPair].first;
    } else if (runsLeft) {
      isPair = false;
      utility = seqHeap.maxUtility();
    } else {
      break;
    }
    
    if(!isPair) {
      Runs best = seqHeap.removeMax();
      if(usedFst[best.symbol] || usedSnd[best.symbol]) continue;
      freqs.decrease(best.symbol, utility + best.frequency);
      if(freqs.getFrequency(currentSymbol) + 1003 >= utility) {
        freqs.increase(best.symbol, utility + best.frequency);
        break;
      }
      usedRun[best.symbol] = true;
      replaceables.push_back(Replacement(
          best.frequency, best.length, best.symbol, false));
      ++currentSymbol;
    } else {
      byte fst = static_cast<byte>((pairs[currentPair].second >> 8) & 0xFF);
      byte snd = static_cast<byte>(pairs[currentPair].second & 0xFF);

      assert(fst != snd);
      if(usedFst[snd] || usedSnd[fst] || usedRun[fst] || usedRun[snd]) {
        ++currentPair;
        continue;
      }
      if(!freqs.decrease(fst, pairs[currentPair].first)) {
        ++currentPair;
        continue;
      }
      if(!freqs.decrease(snd, pairs[currentPair].first)) {
        freqs.increase(fst, pairs[currentPair].first);
        ++currentPair;
        continue;
      }
      
      if(freqs.getFrequency(currentSymbol) + 1003 >= pairs[currentPair].first) {
        freqs.increase(fst, pairs[currentPair].first);
        freqs.increase(snd, pairs[currentPair].first);
        ++currentPair;
        break;
      }
      
      replaceables.push_back(Replacement(
          pairs[currentPair].first, 2, pairs[currentPair].second, true));
      usedFst[fst] = true; usedSnd[snd] = true;
      ++currentPair;
      ++currentSymbol;
    }
  }
}

size_t PairAndRunReplacer::
findEscapeIndex(FrequencyTable& freqs, size_t freeSymbols,
                std::vector<Replacement>& replaceables) {
  if(replaceables.size() <= freeSymbols) return freeSymbols;
  int64 utility = 0;
  size_t i = freeSymbols;
  for(;i < replaceables.size(); ++i) {
    const Replacement& r = replaceables[i];
    if(r.pair) {
      utility += (r.frequency - freqs.getFrequency(i) - 3);
    } else {
      utility += (r.frequency*(r.length-1) - freqs.getFrequency(i) - 3);
    }
  }
  while(utility <= ((int64)freqs.getFrequency(i)) && i > freeSymbols) {
    --i;
    const Replacement& r = replaceables[i];
    if(r.pair) {
      byte fst = static_cast<byte>((r.replaceable & 0xff) >> 8);
      byte snd = static_cast<byte>(r.replaceable & 0xff);
      freqs.increase(fst, r.frequency);
      freqs.increase(snd, r.frequency);
      utility -= (r.frequency - freqs.getFrequency(i) - 3);
    } else {
      int64 score = (r.length - 1)*r.frequency;
      freqs.increase(r.replaceable, score + r.frequency);
      utility -= (score - freqs.getFrequency(i) - 3);
    }
  }
  return i;
}

size_t PairAndRunReplacer::
writeRunReplacement(byte runSymbol, int runLength, byte *dst) const {
  size_t j = 0;
  int pointer = m_runReplacements.listPointer(runSymbol);
  while(pointer != -1 && runLength > 0) {
    const RunReplacement& rr = m_runReplacements.replacement(pointer);
    if((int)rr.length <= runLength) {
      // TODO: optimize (rr.length == 2^k, for some k)
      //size_t times = runLength/rr.length;
      size_t logLen = utils::logFloor(rr.length);
      size_t times = runLength >> logLen;
      std::fill(dst +j, dst + j + times, rr.replacementSymbol);
      runLength -= (times << logLen);
      j += times;
    }
    pointer = rr.nextElement;
  }
  if(runLength > 0) {
    if(m_runReplacements.isEscaped(runSymbol)) {
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


size_t PairAndRunReplacer::
writeReplacedVersion(const byte *src, size_t length, byte* dst) const {
  size_t j = 0;
  byte prev = src[0];
  uint16 pair = prev;
  size_t i = 1;
  while(true) {
    pair = (pair << 8) | src[i]; 
    if(prev != src[i]) {
      byte replValue = m_pairReplacements[pair];
      if(replValue == m_commonByte) {
        dst[j++] = src[i-1];
      } else if(replValue != m_escapeByte) {
        dst[j++] = replValue;
        if(i == length -1) break;
        ++i;
      } else {
        dst[j++] = m_escapeByte;
        dst[j++] = src[i-1];
      }
    } else {
      int runLength = 2;
      ++i;
      while(i < length && prev == src[i]) {++i; ++runLength; }
      j += writeRunReplacement(prev, runLength, dst + j);
    }
    if(i >= length - 1) {
      if(i == length - 1) {
        pair = (src[i] << 8) | 0;
        if(m_pairReplacements[pair] == m_escapeByte && m_escapeByte != m_commonByte)
          dst[j++] = m_escapeByte;
        dst[j++] = src[i];
      }
      break;
    }
    pair = prev = src[i];
    ++i;
  }
  return j;      
}

size_t PairAndRunReplacer::writeHeader(byte* dst) const {
  size_t pos = 0;
  /* First the pairs are written */
  if(m_numOfReplacements == m_runReplacements.size()) {
    dst[0] = dst[1] = dst[2] = 0;
    pos += 3;
  } else {
    byte prevValue = m_escapeByte;
    for(size_t i = 0; i < (1 << 16); ++i) {
      byte val = m_pairReplacements[i];
      if(val != m_commonByte && val != m_escapeByte) {
        prevValue = val;
        dst[pos++] = prevValue;
        dst[pos++] = (i >> 8) & 0xff;
        dst[pos++] = i & 0xff;
      }
    }
    dst[pos++] = prevValue;
    dst[pos++] = (m_escapeByte != m_commonByte)?m_escapeByte:prevValue;
  }
  /* After that runs are written */
  if(m_runReplacements.size() == 0) {
    dst[pos++] = 0; dst[pos++] = 0;
  } else {
      std::vector<std::pair<byte, RunReplacement> > replacements;
      m_runReplacements.replacements(replacements);
      assert(replacements.size() == m_runReplacements.size());
  
      size_t pairs = replacements.size() & 0xfffffffe;
      byte prev = m_escapeByte;
      for(size_t i = 0; i < pairs; i += 2) {
        dst[pos++] = replacements[i].second.replacementSymbol;
        byte lengths = (utils::logFloor((uint32)replacements[i].second.length) << 4) |
            (utils::logFloor((uint32)replacements[i+1].second.length));
        dst[pos++] = lengths;
        dst[pos++] = replacements[i].first;
        dst[pos++] = replacements[i+1].second.replacementSymbol;
        dst[pos++] = replacements[i+1].first;
        prev = replacements[i+1].second.replacementSymbol;
      }

      if(replacements.size() != pairs) {
        prev = replacements.back().second.replacementSymbol;
        dst[pos++] = prev;
        dst[pos++] = utils::logFloor((uint32)replacements.back().second.length) << 4;
        dst[pos++] = replacements.back().first;
        dst[pos++] = (m_runReplacements.escapingInUse())?m_escapeByte:prev;
      } else {
        dst[pos++] = (m_runReplacements.escapingInUse())?m_escapeByte:prev;
        dst[pos++] = 0;
      }
  }
  return pos;
}

size_t PairAndRunReplacer::decideReplacements() {
  FrequencyTable freqTable(m_frequencies);
  size_t freeSymbols = 0;
  while(freqTable.getFrequency(freeSymbols) == 0) ++freeSymbols;

  std::vector<Replacement> replaceables;
  size_t escapeIndex;
  {
    std::vector<std::pair<size_t, uint16> > pairs;
    PairReplacer::makePairList(pairs, m_pairFrequencies);
    if(m_useEscaping) {
      findReplaceablePairsAndRuns(pairs, replaceables, freqTable, 254);
      escapeIndex = (replaceables.size() <= freeSymbols)?freeSymbols:
          findEscapeIndex(freqTable, freeSymbols, replaceables);
    } else { 
      findReplaceablePairsAndRuns(pairs, replaceables, freqTable, freeSymbols);
      escapeIndex = freeSymbols;
    }
  }

  escapeIndex = (replaceables.size() <= freeSymbols)?freeSymbols:
      findEscapeIndex(freqTable, freeSymbols, replaceables);
  m_escapeByte = freqTable.getKey(escapeIndex);

  m_numOfReplacements = (escapeIndex > freeSymbols)?
      escapeIndex:std::min(freeSymbols, replaceables.size());

  if(m_verbose) {
    size_t pairRepls = 0;
    for(size_t i = 0; i < m_numOfReplacements; ++i) {
      if(replaceables[i].pair) ++pairRepls;
    }
    std::clog << "Replacing " << pairRepls << " pairs and "
              << (m_numOfReplacements - pairRepls) << " runs. ";
    if(m_numOfReplacements > freeSymbols) {
      std::clog << "Made " << (escapeIndex - freeSymbols + 1)
                << " symbols free." << std::endl;
    } else {
      std::clog << "No symbols made free." << std::endl;
    }
  }

  m_commonByte = freqTable.getKey(255);
  std::fill(m_pairReplacements, m_pairReplacements + (1 << 16), m_commonByte);

    if(escapeIndex > freeSymbols) {
    for(size_t i = freeSymbols; i <= escapeIndex; ++i) {
      m_runReplacements.addEscaping(freqTable.getKey(i));
      uint16 pairValue = freqTable.getKey(i) << 8;
      for(size_t j = 0; j < 256; ++j, ++pairValue) {
        m_pairReplacements[pairValue] = m_escapeByte;
      }
    }
  }
  if(m_numOfReplacements > 0) {
    for(size_t i = 0; i < m_numOfReplacements; ++i) {
      const Replacement& r = replaceables[i];
      if(r.pair) {
        m_pairReplacements[r.replaceable] = freqTable.getKey(i);
      } else {
        m_runReplacements.addReplacement(r.replaceable, r.length, freqTable.getKey(i));
      }
    }
  } else {
    m_commonByte = m_escapeByte = freqTable.getKey(255);
    std::fill(m_pairReplacements, m_pairReplacements + (1 << 16), m_commonByte);
  }
  m_runReplacements.prepare();

  return m_numOfReplacements;
}

} //namespace pairs_and_runs
} //namespace bwtc
