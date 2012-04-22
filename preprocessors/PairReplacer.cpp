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
    uint32& variables, uint32& specials, uint32& forFree) const
{
  std::sort(pairs.rbegin(), pairs.rend());

  //for(size_t i = 0; i < 256; ++i) {
  //  if(m_grammar.isSpecial(i)) freqs.remove(i);
  //}

  int64 bestUtility = 0;
  FrequencyTable tFreqs(freqs);  
  for(size_t i = 0; i < s_greedyStarts; ++i) {
    FrequencyTable tmpFreqs(tFreqs);
    std::vector<std::pair<uint32, uint16> > tmpReplacements;
    uint32 tVars, tSpecs, tFree;

    int64 utility = findReplaceables(i, pairs, tmpReplacements, tmpFreqs,
                                     maxReplacements, tVars, tSpecs, tFree);

    if(utility > bestUtility) {
      assert(tVars == tmpReplacements.size());
      freqs = tmpFreqs;
      replaceablePairs = tmpReplacements;
      bestUtility = utility;
      variables = tVars;
      specials = tSpecs;
      forFree = tFree;
    }
  }
}

int64 PairReplacer::findReplaceables(size_t startingPair, 
                                     const std::vector<std::pair<uint32, uint16> >& pairs,
                                     std::vector<std::pair<uint32, uint16> >& replPairs,
                                     FrequencyTable& freqs, size_t maxRepl,
                                     uint32& variables, uint32& specials,
                                     uint32& forFree) const {
  bool usedFst[256] = {false}, usedSnd[256] = {true};
  size_t currentPair = startingPair, currentSymbol = 0;
  int64 utility = 0;

  forFree = 0;
  uint32 vars = 0, specs = 0;

  uint32 freeSymbols = 0;
  while(freqs.getFrequency(freeSymbols) == 0) ++freeSymbols;
  
  uint32 withoutNew = m_grammar.specialSymbolPairsLeft() + freeSymbols;
  bool hope = currentPair < pairs.size() && currentSymbol < 256
      && !m_grammar.isSpecial(freqs.getKey(currentSymbol)) && replPairs.size() < maxRepl;

  while(hope) {
    byte fst = static_cast<byte>((pairs[currentPair].second >> 8) & 0xFF);
    byte snd = static_cast<byte>(pairs[currentPair].second & 0xFF);
    
    if(usedFst[snd] || usedSnd[fst] || m_grammar.isSpecial(snd) || m_grammar.isSpecial(fst)) {
      ++currentPair;
      continue;
    }
    /*
    if(m_grammar.isSpecial(freqs.getKey(currentSymbol))) {
      hope = false;
      break;
      }*/
    
    freqs.decrease(fst, pairs[currentPair].first);
    freqs.decrease(snd, pairs[currentPair].first);

    uint32 fr = freqs.getFrequency(currentSymbol);

    if(freqs.getFrequencyWithKey(fst) == 0) {
      ++withoutNew;
    }
    if(fst != snd && freqs.getFrequencyWithKey(snd) == 0) {
      ++withoutNew;
    }


    if(fr + 1003 >= pairs[currentPair].first || withoutNew == 0) {
      freqs.increase(fst, pairs[currentPair].first);
      freqs.increase(snd, pairs[currentPair].first);
      hope = false;
      break;
    }

    --withoutNew;
    if(fr == 0) ++forFree;

    utility += pairs[currentPair].first;
    utility -= fr;

    replPairs.push_back(pairs[currentPair]);
    ++vars;
    ++currentPair;
    ++currentSymbol;

    usedFst[fst] = true;
    usedSnd[snd] = true;
    hope = currentPair < pairs.size() && currentSymbol < 256
        && !m_grammar.isSpecial(freqs.getKey(currentSymbol));
    
  }

  hope = currentPair < pairs.size() && currentSymbol < 256
      && !m_grammar.isSpecial(freqs.getKey(currentSymbol)) && withoutNew == 0;

  uint32 tSpecials = m_grammar.numberOfSpecialSymbols();
  uint32 unusedSpecials = 0;

  while(hope) {
    uint32 newFrees = 0;
    int64 utilityAfterNewSpecials = 0;
    int64 utilityFromFrees = 0;
    int64 utilityFromFreesBeg = 0;
    //Take new special
    uint32 limit = 2*tSpecials + unusedSpecials;
    uint32 nSpecials = 1;
    //while(currentSymbol < 255 && m_grammar.isSpecial(freqs.getKey(currentSymbol))) ++currentSymbol;
    utilityAfterNewSpecials -= freqs.getFrequency(currentSymbol);
    std::vector<byte> sSymbols;
    sSymbols.push_back(freqs.getKey(currentSymbol));
    assert(freqs.getFrequency(currentSymbol) > 0);
    //freqs.removeKey(currentSymbol++);
    currentSymbol++;
    if(tSpecials == 0) {
      limit = 2;
      nSpecials = 2;
      //while(currentSymbol < 255 && m_grammar.isSpecial(freqs.getKey(currentSymbol))) ++currentSymbol;
      utilityAfterNewSpecials -= freqs.getFrequency(currentSymbol);
      sSymbols.push_back(freqs.getKey(currentSymbol));
      //freqs.removeKey(currentSymbol++);
      currentSymbol++;
    }

    uint32 pairsAfterSpecial = 0;
    uint32 freeVars = 0;
    uint32 freeVarsBeg = 0;
    std::vector<std::pair<byte, byte> > tPairs;
    while(pairsAfterSpecial < limit && hope) {
      /*
        if(m_grammar.isSpecial(freqs.getKey(currentSymbol))) {
        hope = false;
        break;
        }*/

      byte fst = static_cast<byte>((pairs[currentPair].second >> 8) & 0xFF);
      byte snd = static_cast<byte>(pairs[currentPair].second & 0xFF);

      if(usedFst[snd] || usedSnd[fst] || m_grammar.isSpecial(snd) || m_grammar.isSpecial(fst)) {
        ++currentPair;
        continue;
      }

      freqs.decrease(fst, pairs[currentPair].first);
      freqs.decrease(snd, pairs[currentPair].first);

      if(freqs.getFrequencyWithKey(fst) == 0) {
        ++newFrees;
      }
      if(fst != snd && freqs.getFrequencyWithKey(snd) == 0) ++newFrees;

      uint32 fr = freqs.getFrequency(currentSymbol);
      if(newFrees > 0) fr = 0;

      
      if(fr + 1003 >= pairs[currentPair].first) {
        freqs.increase(fst, pairs[currentPair].first);
        freqs.increase(snd, pairs[currentPair].first);
        hope = false;
        break;
      }

      bool reallyFree = true;
      if(fst == sSymbols[0] || snd == sSymbols[0]) reallyFree = false;
      if(sSymbols.size() > 1 && (fst == sSymbols[1] || snd == sSymbols[1]))
         reallyFree = false;

      if(newFrees > 0) {
        ++forFree;
        --newFrees;
        if(reallyFree && pairsAfterSpecial == 0) {
          ++freeVarsBeg;
          utilityFromFreesBeg += pairs[currentPair].first;
        } else {
          ++freeVars;
          utilityFromFrees += pairs[currentPair].first;
        }
      } else {
        utilityAfterNewSpecials += pairs[currentPair].first;
        utilityAfterNewSpecials -= fr;
        ++pairsAfterSpecial;
      }
      
      replPairs.push_back(pairs[currentPair]);
      ++currentPair;
      ++currentSymbol;


      usedFst[fst] = true;
      usedSnd[snd] = true;

      hope = currentPair < pairs.size() && currentSymbol < 256 &&
          !m_grammar.isSpecial(freqs.getKey(currentSymbol));
    }
    // Worth of new special symbol
    int64 totalUt = utilityFromFrees + utilityAfterNewSpecials + utilityFromFreesBeg;
    

    if(totalUt > (freeVars+freeVarsBeg+pairsAfterSpecial)*1000 && pairsAfterSpecial > freeVarsBeg) {
      utility += totalUt;
      vars += pairsAfterSpecial + freeVars + freeVarsBeg;
      specs += nSpecials;
      tSpecials += nSpecials;
      unusedSpecials = freeVars + freeVarsBeg;
    } else if(utilityFromFreesBeg > 0){
      utility += utilityFromFreesBeg;
      vars += freeVarsBeg;
      forFree -= freeVars;
      hope = false;
      for(size_t i = 0; i < pairsAfterSpecial + freeVars; ++i) {
        byte fst = static_cast<byte>((replPairs.back().second >> 8) & 0xFF);
        byte snd = static_cast<byte>(replPairs.back().second & 0xFF);

        freqs.increase(fst, replPairs.back().first);
        freqs.increase(snd, replPairs.back().first);
        assert(!replPairs.empty());
        replPairs.pop_back();
      }
    } else {
      forFree = forFree - freeVarsBeg - freeVars;
      hope = false;
      for(size_t i = 0; i < pairsAfterSpecial + freeVars + freeVarsBeg; ++i) {
        byte fst = static_cast<byte>((replPairs.back().second >> 8) & 0xFF);
        byte snd = static_cast<byte>(replPairs.back().second & 0xFF);

        freqs.increase(fst, replPairs.back().first);
        freqs.increase(snd, replPairs.back().first);
        assert(!replPairs.empty());
        replPairs.pop_back();
      }
    }
  }
  specials = specs;
  variables = vars;
  return utility;
}

void PairReplacer::constructReplacementTable(
    const std::vector<std::pair<uint32, uint16> >& pairs,
    const std::vector<byte>& freedSymbols,
    const std::vector<byte>& newSpecials,
    const std::vector<byte>& replacements)
{
  m_grammar.beginUpdatingRules();
  assert(m_numOfReplacements > 0);
  assert(m_numOfReplacements == replacements.size());
  for(size_t i = 0; i < m_numOfReplacements; ++i) {
    assert(m_commonByte != replacements[i]);
    m_replacements[pairs[i].second] = (replacements[i]<<8) | m_commonByte;
    m_grammar.addRule(replacements[i], (pairs[i].second >> 8) & 0xff,
                      pairs[i].second & 0xff);
  }
  std::vector<uint16> nextSpecialPairs;
  m_grammar.expandAlphabet(freedSymbols, newSpecials, nextSpecialPairs);
  //replacements.size()- freedSymbols.size());
  assert(freedSymbols.size() == nextSpecialPairs.size());

  for(size_t i = 0; i < freedSymbols.size(); ++i) {
    assert(m_commonByte != freedSymbols[i]);
    uint16 special = nextSpecialPairs[i];
    uint16 hVal = freedSymbols[i] << 8;
    for(size_t j = 0; j < 256; ++j) {
      if((m_replacements[hVal | j] >> 8) == m_commonByte)
        m_replacements[hVal | j] = special;
    }
  }
  for(size_t i = 0; i < newSpecials.size(); ++i) {
    assert(m_commonByte != newSpecials[i]);
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

  makePairList(pairs, m_pairFrequencies);
  std::vector<std::pair<uint32, uint16> > replaceablePairs;

  uint32 nSpecials = 0;
  uint32 nVariables = 0;
  uint32 forFree = 0;
  
  findReplaceablePairs(pairs, replaceablePairs, freqTable, 254,
                       nVariables, nSpecials, forFree);


  assert(nVariables == replaceablePairs.size());
  m_numOfNewSpecials = nSpecials;
  m_numOfFreedSymbols = (forFree > nVariables)?0:(nVariables - forFree);

  std::vector<byte> newSpecials;
  std::vector<byte> freedSymbols;
  std::vector<byte> replacements;

  assert(nVariables + m_numOfNewSpecials < 256);

  size_t j = 0;
  // 'Free' variables must be gathered first, but they need to be
  // used last for replacements for the correct interpretation
  // of rules
  std::vector<byte> frees;
  for(size_t i = 0; i < std::min(forFree,nVariables); ++i) {
    byte k = freqTable.getKey(j);
    assert(!m_grammar.isSpecial(k));
    assert(freqTable.getFrequency(j) == 0);
    frees.push_back(k);
    ++j;
  }
  for(size_t i = 0; i < m_numOfNewSpecials; ++i) {
    byte k = freqTable.getKey(j++);
    assert(!m_grammar.isSpecial(k));
    newSpecials.push_back(k);
  }
  for(size_t i = 0; i < m_numOfFreedSymbols; ++i) {
    byte k = freqTable.getKey(j++);
    assert(!m_grammar.isSpecial(k));
    freedSymbols.push_back(k);
    replacements.push_back(k);
  }
  for(size_t i = 0; i < frees.size(); ++i) {
    replacements.push_back(frees[i]);
  }

  assert(replacements.size() == replaceablePairs.size());
  
  m_numOfReplacements = replaceablePairs.size();

  if(m_verbose) {
    std::clog << "Replacing " << m_numOfReplacements << " pairs. ";
    std::clog << "Freed " << m_numOfFreedSymbols << " symbols for pairs and "
              << m_numOfNewSpecials << " symbols for special "
              << "characters." << std::endl;
  }
  
  m_commonByte = freqTable.getKey(255);
  std::fill(m_replacements, m_replacements + (1 << 16), (m_commonByte<<8)|m_commonByte);
  if(m_numOfReplacements > 0)
    constructReplacementTable(replaceablePairs,freedSymbols, newSpecials,replacements);
  return m_numOfReplacements;
}

} //namespace bwtc
