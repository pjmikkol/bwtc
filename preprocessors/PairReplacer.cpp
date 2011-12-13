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

PairReplacer::PairReplacer()
    : m_prev(0), m_numOfReplacements(0), m_escapeByte(0), m_commonByte(0),
      m_analysationStarted(false), m_verbose(false)
{}

PairReplacer::PairReplacer(bool verbose)
    : m_prev(0), m_numOfReplacements(0), m_escapeByte(0), m_commonByte(0),
      m_analysationStarted(false), m_verbose(verbose)      
{}

PairReplacer::PairReplacer(const PairReplacer& pr)
    : m_prev(pr.m_prev), m_numOfReplacements(pr.m_numOfReplacements),
      m_escapeByte(pr.m_escapeByte),
      m_analysationStarted(pr.m_analysationStarted),
      m_verbose(pr.m_verbose)
{
  std::copy(pr.m_frequencies, pr.m_frequencies + 256, m_frequencies);
  std::copy(pr.m_replacements, pr.m_replacements + (1 << 16), m_replacements);
  std::copy(pr.m_pairFrequencies, pr.m_pairFrequencies + (1 << 16),
            m_pairFrequencies);
}

PairReplacer::~PairReplacer() {}

void PairReplacer::analyseData(const byte *data, size_t length, bool reset) {
  assert(length > 0);
  size_t i = 0;
  if(!m_analysationStarted) beginAnalysing(data[i++], reset);

  for(; i < length; ++i) {
    m_prev = (m_prev << 8) | data[i];

    ++m_pairFrequencies[m_prev];
    ++m_frequencies[data[i]];
  }
}

void PairReplacer::finishAnalysation() {}

void PairReplacer::analyseData(byte next) {
  assert(m_analysationStarted);
  m_prev = (m_prev << 8) | next;
  
  ++m_pairFrequencies[m_prev];
  ++m_frequencies[next];
}

void PairReplacer::beginAnalysing(byte first, bool reset) {
  assert(!m_analysationStarted);
  if(reset) {
    std::fill(m_pairFrequencies, m_pairFrequencies + (1 << 16), 0);
    std::fill(m_frequencies, m_frequencies + 256, 0);
  }

  m_analysationStarted = true;
  m_prev = first;
  ++m_frequencies[first];
}

void PairReplacer::makePairList(std::vector<std::pair<size_t, uint16> >& pairs,
                                const size_t *pairFrequencies) const
{
  assert(pairs.empty());
  for(size_t i = 0; i < (1 << 16); ++i) {
    if((i & 0xff) == ((i >> 8) & 0xff)) continue;
    pairs.push_back(std::make_pair(pairFrequencies[i], i));
  }
}

void PairReplacer::findReplaceablePairs(
    std::vector<std::pair<size_t, uint16> >& pairs,
    std::vector<std::pair<size_t, uint16> >& replaceablePairs,
    FrequencyTable& freqs) const
{
  bool usedFst[256] = {false}, usedSnd[256] = {true};
  size_t currentPair = 0, currentSymbol = 0;

  std::sort(pairs.rbegin(), pairs.rend());

  while(replaceablePairs.size() < 254 && currentPair < pairs.size() &&
        currentSymbol < 254)
  {
    byte fst = static_cast<byte>((pairs[currentPair].second >> 8) & 0xFF);
    byte snd = static_cast<byte>(pairs[currentPair].second & 0xFF);

    assert(fst != snd);
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

    if(freqs.getFrequency(currentSymbol) + 3 >= pairs[currentPair].first) {
      freqs.increase(fst, pairs[currentPair].first);
      freqs.increase(snd, pairs[currentPair].first);
      ++currentPair;
      break;
    }

    replaceablePairs.push_back(pairs[currentPair]);
    usedFst[fst] = true; usedSnd[snd] = true;
    ++currentPair;
    ++currentSymbol;
  }
  assert(currentSymbol == replaceablePairs.size());
}

size_t PairReplacer::findEscapeIndex(FrequencyTable& freqs, size_t freeSymbols,
                                     std::vector<std::pair<size_t, uint16> >&
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
    const std::vector<std::pair<size_t, uint16> >& pairs,
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

size_t PairReplacer::writeHeader(byte *to) const {
  if(m_numOfReplacements == 0)  {
    to[0] = to[1] = to[2] = 0;
    return 3;
  }
  size_t pos = 0;
  byte prevValue = m_escapeByte;
  for(size_t i = 0; i < (1 << 16); ++i) {
    byte val = m_replacements[i];
    if(val != m_commonByte && val != m_escapeByte) {
      prevValue = val;
      to[pos++] = prevValue;
      to[pos++] = (i >> 8) & 0xff;
      to[pos++] = i & 0xff;
    }
  }
  to[pos++] = prevValue;
  to[pos++] = (m_escapeByte != m_commonByte)?m_escapeByte:prevValue;
  return pos;    
}

size_t PairReplacer::writeReplacedVersion(const byte *src, size_t length, byte *dst) const
{
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
  std::vector<std::pair<size_t, uint16> > pairs;

  size_t freeSymbols = 0;
  while(freqTable.getFrequency(freeSymbols) == 0) ++freeSymbols;

  makePairList(pairs, m_pairFrequencies);
  std::vector<std::pair<size_t, uint16> > replaceablePairs;
  findReplaceablePairs(pairs, replaceablePairs, freqTable);


  size_t escapeIndex = (replaceablePairs.size() > freeSymbols)?
      findEscapeIndex(freqTable, freeSymbols, replaceablePairs):freeSymbols;

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
  
  if(m_numOfReplacements > 0) {
    m_commonByte = freqTable.getKey(255);
    std::fill(m_replacements, m_replacements + (1 << 16), m_commonByte);
    m_escapeByte =(m_numOfReplacements <= freeSymbols)?
        m_commonByte:freqTable.getKey(escapeIndex);
    constructReplacementTable(replaceablePairs, freqTable, freeSymbols);
  } else {
    m_commonByte = m_escapeByte = freqTable.getKey(255);
    std::fill(m_replacements, m_replacements + (1 << 16), m_commonByte);
  }

  return m_numOfReplacements;
}

} //namespace bwtc
