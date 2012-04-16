/**
 * @file PairReplacer.cpp
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
 * Implementation for preprocessor which replaces the most often
 * occurring pairs.
 */
#include "PairReplacer.hpp"
#include "FrequencyTable.hpp"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <vector>

namespace bwtc {

PairReplacer::PairReplacer(Grammar& grammar, bool useEscaping)
    : m_grammar(grammar), m_prev(0), m_numOfReplacements(0), m_escapeByte(0),
      m_commonByte(0), m_analysationStarted(false), m_verbose(false),
      m_useEscaping(useEscaping)
{}

PairReplacer::PairReplacer(Grammar& grammar, bool useEscaping, bool verbose)
    : m_grammar(grammar), m_prev(0), m_numOfReplacements(0), m_escapeByte(0),
      m_commonByte(0), m_analysationStarted(false), m_verbose(verbose),
      m_useEscaping(useEscaping)
{}

PairReplacer::~PairReplacer() {}

void PairReplacer::analyseData(const byte *data, size_t length, bool reset) {
  assert(length > 0);
  size_t i = 0;
  if(!m_analysationStarted) beginAnalysing(data[i++], reset);

  /*
  for(; i < length; ++i) {
    analyseData(data[i]);
    }*/
  for(; i < (length & 0xfffffffe); ++i) {
    analyseData0(data[i]);
    analyseData(data[++i]);
  }
  if(length & 0xfffffffe) analyseData0(data[i]);
}

void PairReplacer::finishAnalysation() {}

void PairReplacer::resetAnalyseData() {
  std::fill(m_pairFrequencies, m_pairFrequencies + (1 << 16), 0);
  std::fill(m_frequencies, m_frequencies + 256, 0);
}

void PairReplacer::beginAnalysing(byte first, bool reset) {
  beginAnalysing(reset);
  m_prev = first;
  ++m_frequencies[first];
}

void PairReplacer::beginAnalysing(bool reset) {
  assert(!m_analysationStarted);
  if(reset) resetAnalyseData();

  m_analysationStarted = true;
}

void PairReplacer::
makePairList(std::vector<std::pair<uint32, uint16> >& pairs,
             const size_t *pairFrequencies) {
  assert(pairs.empty());
  for(size_t i = 0; i < (1 << 16); ++i) {
    //if((i & 0xff) == ((i >> 8) & 0xff)) continue;
    pairs.push_back(std::make_pair(pairFrequencies[i], i));
  }
}

void PairReplacer::findReplaceablePairs(
    std::vector<std::pair<uint32, uint16> >& pairs,
    std::vector<std::pair<uint32, uint16> >& replaceablePairs,
    FrequencyTable& freqs, size_t maxReplacements) const
{
  std::sort(pairs.rbegin(), pairs.rend());

  int64 bestUtility = 0;
  FrequencyTable bestFreqs(freqs);
  std::vector<std::pair<uint32, uint16> > bestReplacements;
  for(size_t i = 0; i < s_greedyStarts; ++i) {
    FrequencyTable tmpFreqs(freqs);
    std::vector<std::pair<uint32, uint16> > tmpReplacements;
    int64 utility = findReplaceablePairs(i, pairs, tmpReplacements, tmpFreqs, maxReplacements);
    if(utility > bestUtility) {
      bestFreqs = tmpFreqs;
      bestReplacements = tmpReplacements;
      bestUtility = utility;
    }
  }
  freqs = bestFreqs;
  replaceablePairs = bestReplacements;
}

int64 PairReplacer::findReplaceablePairs(
    size_t startingPair,
    const std::vector<std::pair<uint32, uint16> >& pairs,
    std::vector<std::pair<uint32, uint16> >& replaceablePairs,
    FrequencyTable& freqs, size_t maxReplacements) const
{
  bool usedFst[256] = {false}, usedSnd[256] = {true};
  size_t currentPair = startingPair, currentSymbol = 0;

  int64 utility = 0;

  while(replaceablePairs.size() < maxReplacements && currentPair < pairs.size()
        && currentSymbol < maxReplacements)
  {
    byte fst = static_cast<byte>((pairs[currentPair].second >> 8) & 0xFF);
    byte snd = static_cast<byte>(pairs[currentPair].second & 0xFF);

    if(usedFst[snd] || usedSnd[fst]) {
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

    uint32 fr = freqs.getFrequency(currentSymbol);
    if(fr + 1003 >= pairs[currentPair].first) {
      freqs.increase(fst, pairs[currentPair].first);
      freqs.increase(snd, pairs[currentPair].first);
      ++currentPair;
      break;
    }
    utility += pairs[currentPair].first;
    utility -= fr;

    replaceablePairs.push_back(pairs[currentPair]);
    usedFst[fst] = true; usedSnd[snd] = true;
    ++currentPair;
    ++currentSymbol;
  }
  assert(currentSymbol == replaceablePairs.size());
  return utility;
}

size_t PairReplacer::findEscapeIndex(FrequencyTable& freqs, size_t freeSymbols,
                                     std::vector<std::pair<uint32, uint16> >&
                                     suitablePairs)
{
  if(suitablePairs.size() <= freeSymbols) return freeSymbols;
  int64 utility = 0;
  size_t i = freeSymbols;
  for(;i < suitablePairs.size(); ++i) {
    utility += (suitablePairs[i].first - freqs.getFrequency(i) - 3);
  }
  while(utility <= static_cast<int64>(freqs.getFrequency(i)) && i > freeSymbols)
  {
    --i;
    byte fst = static_cast<byte>((suitablePairs[i].second & 0xff) >> 8);
    byte snd = static_cast<byte>(suitablePairs[i].second & 0xff);
    freqs.increase(fst, suitablePairs[i].first);
    freqs.increase(snd, suitablePairs[i].first);
    utility -= (suitablePairs[i].first - freqs.getFrequency(i) - 3);
  }
  return i;
}

void PairReplacer::constructReplacementTable(
    const std::vector<std::pair<uint32, uint16> >& pairs,
    const FrequencyTable& freqTable, size_t freeSymbols)
{
  assert(m_numOfReplacements > 0);
  if(m_escapeByte != m_commonByte) {
    for(int i = freeSymbols; i < m_numOfReplacements+1; ++i) {
      uint16 pairValue = freqTable.getKey(i) << 8;
      for(size_t j = 0; j < 256; ++j, ++pairValue) {
        m_replacements[pairValue] = m_escapeByte;
      }
    }
  }
  for(size_t i = 0; i < m_numOfReplacements; ++i) {
    m_replacements[pairs[i].second] = freqTable.getKey(i);
  }
}

size_t PairReplacer::
writeReplacedVersion(const byte *src, size_t length, byte *dst) const {
  size_t j = 0; /* Is used for indexing the target. */
  uint16 pair = src[0];
  size_t i = 1;
  while(true) {
    pair = (pair << 8) | src[i];
    byte replValue = m_replacements[pair];
    if(replValue == m_commonByte) {
      dst[j++] = src[i-1];
    } else if(replValue != m_escapeByte) {
      dst[j++] = replValue;
      if(i == length - 1) break;
      pair = src[++i];
    } else {
      dst[j++] = m_escapeByte;
      dst[j++] = src[i-1];
    }
    if(i == length - 1) {
      pair = (src[i] << 8) | 0;
      if(m_replacements[pair] == m_escapeByte && m_escapeByte != m_commonByte)
        dst[j++] = m_escapeByte;
      dst[j++] = src[i];
      break;
    }
    ++i;
  }
  return j;
}

size_t PairReplacer::decideReplacements() {
  assert(m_analysationStarted);
  FrequencyTable freqTable(m_frequencies);
  std::vector<std::pair<uint32, uint16> > pairs;

  size_t freeSymbols = 0;
  while(freqTable.getFrequency(freeSymbols) == 0) ++freeSymbols;

  makePairList(pairs, m_pairFrequencies);
  std::vector<std::pair<uint32, uint16> > replaceablePairs;

  size_t escapeIndex;
  
  if(m_useEscaping) {
    findReplaceablePairs(pairs, replaceablePairs, freqTable, 254);
    escapeIndex = (replaceablePairs.size() > freeSymbols)?
        findEscapeIndex(freqTable, freeSymbols, replaceablePairs):freeSymbols;
  } else {
    if(freeSymbols > 0)
      findReplaceablePairs(pairs, replaceablePairs, freqTable, freeSymbols);
    escapeIndex = freeSymbols;
  }

  m_numOfReplacements = (escapeIndex > freeSymbols)?
      escapeIndex: std::min(freeSymbols, replaceablePairs.size());
  if(m_verbose) {
    std::clog << "Replacing " << m_numOfReplacements << " pairs. ";
    if(m_numOfReplacements > freeSymbols)
      std::clog << "Made " << (escapeIndex - freeSymbols + 1)
                << " symbols free." << std::endl;
    else
      std::clog << "No symbols made free." << std::endl;
  }
  
  m_commonByte = freqTable.getKey(255);
  std::fill(m_replacements, m_replacements + (1 << 16), m_commonByte);
  if(m_numOfReplacements > 0) {
    m_escapeByte = (m_numOfReplacements <= freeSymbols)?
        m_commonByte:freqTable.getKey(escapeIndex);
    constructReplacementTable(replaceablePairs, freqTable, freeSymbols);
  } else {
    m_escapeByte = m_commonByte;
  }

  return m_numOfReplacements;
}

} //namespace bwtc
