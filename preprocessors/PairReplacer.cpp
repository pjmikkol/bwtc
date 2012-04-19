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
#include "Grammar.hpp"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <vector>

namespace bwtc {

PairReplacer::PairReplacer(Grammar& grammar)
    : m_grammar(grammar), m_prev(0), m_numOfReplacements(0), m_numOfFreedSymbols(0),
      m_numOfNewSpecials(0), m_commonByte(0), m_analysationStarted(false),
      m_verbose(false)
{}

PairReplacer::PairReplacer(Grammar& grammar, bool verbose)
    : m_grammar(grammar), m_prev(0), m_numOfReplacements(0), m_numOfFreedSymbols(0),
      m_numOfNewSpecials(0), m_commonByte(0), m_analysationStarted(false),
      m_verbose(verbose)
{}

PairReplacer::~PairReplacer() {}

void PairReplacer::analyseData(const byte *data, size_t length, bool reset) {
  assert(length > 2);
  size_t i = 0;
  if(!m_analysationStarted) beginAnalysing(data[i++], reset);

  /*
  for(; i < length; ++i) {
    analyseData(data[i]);
    }*/
  for(; i < ((length-1) & 0xfffffffe); ++i) {
    analyseData0(data[i]);
    analyseData(data[++i]);
  }
  if((length & 0x1) == 0) analyseData0(data[i]);
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
  const uint32* grammarFreq = m_grammar.frequencies();
  for(int i = 0; i < 256; ++i) m_frequencies[i] += grammarFreq[i];
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
    pairs.push_back(std::make_pair(pairFrequencies[i], i));
  }
}

void PairReplacer::findReplaceablePairs(
    std::vector<std::pair<uint32, uint16> >& pairs,
    std::vector<std::pair<uint32, uint16> >& replaceablePairs,
    FrequencyTable& freqs, size_t maxReplacements,
    std::vector<byte>& freedSymbols,
    std::vector<byte>& newSpecials,
    std::vector<byte>& replacements) const
{
  std::sort(pairs.rbegin(), pairs.rend());

  int64 bestUtility = 0;
  for(size_t i = 0; i < s_greedyStarts; ++i) {
    FrequencyTable tmpFreqs(freqs);
    std::vector<std::pair<uint32, uint16> > tmpReplacements;
    std::vector<byte> tmpFreedSymbols;
    std::vector<byte> tmpNewSpecials;
    std::vector<byte> tmpVariables;
    int64 utility = findReplaceablePairs(i, pairs, tmpReplacements, tmpFreqs,
                                         maxReplacements, tmpFreedSymbols,
                                         tmpNewSpecials, tmpVariables);
    if(utility > bestUtility) {
      freqs = tmpFreqs;
      replaceablePairs = tmpReplacements;
      bestUtility = utility;
      newSpecials = tmpNewSpecials;
      freedSymbols = tmpFreedSymbols;
      replacements = tmpVariables;
    }
  }
}

int64 PairReplacer::findReplaceablePairs(
    size_t startingPair,
    const std::vector<std::pair<uint32, uint16> >& pairs,
    std::vector<std::pair<uint32, uint16> >& replaceablePairs,
    FrequencyTable& freqs, size_t maxReplacements,
    std::vector<byte>& freedSymbols,
    std::vector<byte>& newSpecials,
    std::vector<byte>& replacements) const
{
  bool usedFst[256] = {false}, usedSnd[256] = {true};
  size_t currentPair = startingPair, currentSymbol = 0;

  int64 utility = 0;

  uint32 freeSpecialPairs = m_grammar.specialSymbolPairsLeft();
  bool noHope = false;

  while(replaceablePairs.size() < maxReplacements && currentPair < pairs.size()
        && currentSymbol < maxReplacements)
  {
    byte fst = static_cast<byte>((pairs[currentPair].second >> 8) & 0xFF);
    byte snd = static_cast<byte>(pairs[currentPair].second & 0xFF);
    
    if(usedFst[snd] || usedSnd[fst] || m_grammar.isSpecial(snd) || m_grammar.isSpecial(fst)) {
      ++currentPair;
      continue;
    }
    if(m_grammar.isSpecial(freqs.getKey(currentSymbol))) {
      ++currentSymbol;
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
    if(fr > 0 && freeSpecialPairs == 0) {
      freqs.increase(fst, pairs[currentPair].first);
      freqs.increase(snd, pairs[currentPair].first);
      break;
    }

    if(fr + 1003 >= pairs[currentPair].first) {
      freqs.increase(fst, pairs[currentPair].first);
      freqs.increase(snd, pairs[currentPair].first);
      ++currentPair;
      noHope = true;
      break;
    }
    if(fr > 0) {
      --freeSpecialPairs;
      freedSymbols.push_back(freqs.getKey(currentSymbol));
    }
    
    utility += pairs[currentPair].first;
    utility -= fr;
    

    replaceablePairs.push_back(pairs[currentPair]);
    usedFst[fst] = true; usedSnd[snd] = true;
    replacements.push_back(freqs.getKey(currentSymbol));
    ++currentPair;
    ++currentSymbol;
  }

  assert(currentSymbol == replaceablePairs.size());
  if(noHope) return utility;
  uint32 specialSymbols = m_grammar.numberOfSpecialSymbols();

  int64 totalCostAfterNewSpecial = 0;
  while(replaceablePairs.size() < maxReplacements && currentPair < pairs.size()
        && currentSymbol < maxReplacements)
  {
    while(m_grammar.isSpecial(currentSymbol) && currentSymbol < maxReplacements) {
      ++currentSymbol;
    }
    newSpecials.push_back(freqs.getKey(currentSymbol));
    //Make currentSymbol to new special symbol
    totalCostAfterNewSpecial -= freqs.getFrequency(currentSymbol);
    // one special pair goes for the original occurrences of the new special
    // character
    freeSpecialPairs = 2*specialSymbols++;
    ++currentSymbol;

    std::vector<byte> tmpFreedSymbols;
    std::vector<byte> tmpReplacements;
    std::vector<std::pair<uint32, uint16> > tmpReplaceablePairs;
    while(replaceablePairs.size() + tmpReplaceablePairs.size() < maxReplacements
          && currentPair < pairs.size()
          && currentSymbol < maxReplacements && freeSpecialPairs > 0) {
      byte fst = static_cast<byte>((pairs[currentPair].second >> 8) & 0xFF);
      byte snd = static_cast<byte>(pairs[currentPair].second & 0xFF);
      
      if(usedFst[snd] || usedSnd[fst] || m_grammar.isSpecial(snd) || m_grammar.isSpecial(fst)) {
        ++currentPair;
        continue;
      }
      if(m_grammar.isSpecial(freqs.getKey(currentSymbol))) {
        ++currentSymbol;
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
      
      totalCostAfterNewSpecial += pairs[currentPair].first;
      totalCostAfterNewSpecial -= fr;
      --freeSpecialPairs;
      tmpFreedSymbols.push_back(freqs.getKey(currentSymbol));
      tmpReplacements.push_back(freqs.getKey(currentSymbol));

      tmpReplaceablePairs.push_back(pairs[currentPair]);
      usedFst[fst] = true; usedSnd[snd] = true;
      ++currentPair;
      ++currentSymbol;
    }
    
    if(specialSymbols > 1 &&
       totalCostAfterNewSpecial < (int64)tmpReplaceablePairs.size()*1000) {
      newSpecials.pop_back();
      if(newSpecials.size() == 1 && m_grammar.numberOfSpecialSymbols() == 0)
        newSpecials.pop_back();
      break;
    } else {
      for(size_t i = 0; i < tmpFreedSymbols.size(); ++i) {
        freedSymbols.push_back(tmpFreedSymbols[i]);
        replaceablePairs.push_back(tmpReplaceablePairs[i]);
        replacements.push_back(tmpReplacements[i]);
      }
      utility += totalCostAfterNewSpecial;
      totalCostAfterNewSpecial = 0;
    }
    
  }
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
    const std::vector<byte>& freedSymbols,
    const std::vector<byte>& newSpecials,const std::vector<byte>& replacements)
{
  m_grammar.beginUpdatingRules();
  assert(m_numOfReplacements > 0);
  assert(m_numOfReplacements == replacements.size());
  for(size_t i = 0; i < m_numOfReplacements; ++i) {
    m_replacements[pairs[i].second] = (replacements[i]<<8) | m_commonByte;
    m_grammar.addRule(replacements[i], (pairs[i].second >> 8) & 0xff, pairs[i].second & 0xff);
  }
  std::vector<uint16> nextSpecialPairs;
  m_grammar.expandAlphabet(freedSymbols, newSpecials, nextSpecialPairs);

  for(size_t i = 0; i < freedSymbols.size(); ++i) {
    uint16 special = nextSpecialPairs[i];
    uint16 hVal = freedSymbols[i] << 8;
    for(size_t j = 0; j < 256; ++j) {
      if((m_replacements[hVal | j] >> 8) == m_commonByte)
        m_replacements[hVal | j] = special;
    }
  }
  for(size_t i = 0; i < newSpecials.size(); ++i) {
    uint16 special = (newSpecials[i] << 8) | newSpecials[i];
    uint16 hVal = newSpecials[i] << 8;
    for(size_t j = 0; j < 256; ++j) {
      if((m_replacements[hVal | j] >> 8) == m_commonByte)
        m_replacements[hVal | j] = special;
    }
  }
  m_grammar.endUpdatingRules();
}

size_t PairReplacer::
writeReplacedVersion(const byte *src, size_t length, byte *dst) const {
  size_t j = 0; /* Is used for indexing the target. */
  uint16 pair = src[0];
  size_t i = 1;
  uint16 noop = (m_commonByte << 8) | m_commonByte;
  while(true) {
    pair = (pair << 8) | src[i];
    uint16 replValue = m_replacements[pair];
    if(replValue == noop) {
      dst[j++] = src[i-1];
    } else if((replValue & 0xff) == m_commonByte) {
      dst[j++] = replValue >> 8;
      if(i == length - 1) break;
      pair = src[++i];
    } else {
      dst[j++] = replValue >> 8;
      dst[j++] = replValue & 0xff;
    }
    if(i == length - 1) {
      pair = (src[i] << 8) | 0;
      if(((replValue = m_replacements[pair]) & 0xff) != m_commonByte) {
        dst[j++] = replValue >> 8;
        dst[j++] = replValue & 0xff;
      } else {
        dst[j++] = src[i];
      }
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

  std::vector<byte> replacements;
  std::vector<byte> freedSymbols;
  std::vector<byte> newSpecials;

  findReplaceablePairs(pairs, replaceablePairs, freqTable, 254,
                       freedSymbols, newSpecials, replacements);
  m_numOfFreedSymbols = freedSymbols.size();
  m_numOfNewSpecials = newSpecials.size();

  m_numOfReplacements = replaceablePairs.size();

  if(m_verbose) {
    std::clog << "Replacing " << m_numOfReplacements << " pairs. ";
    std::clog << "Freed " << m_numOfFreedSymbols << " symbols for pairs and "
              << m_numOfNewSpecials << " symbols for special "
              << "characters." << std::endl;
  }
  
  m_commonByte = freqTable.getKey(255);
  std::fill(m_replacements, m_replacements + (1 << 16), (m_commonByte<<8)|m_commonByte);
  constructReplacementTable(replaceablePairs,freedSymbols, newSpecials,replacements);
  return m_numOfReplacements;
}

} //namespace bwtc
