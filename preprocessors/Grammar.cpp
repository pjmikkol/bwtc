/**
 * @file Grammar.cpp
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
 * Implementation for grammar-object used during the preprocessing.
 */

#include <algorithm>

#include "../globaldefs.hpp"
#include "../Utils.hpp"
#include "Grammar.hpp"

namespace bwtc {

Grammar::Grammar() : m_specialSymbolsAsVariables(0), m_updatingRules(false) {
  std::fill(m_isSpecialSymbol, m_isSpecialSymbol + 256, false);
  std::fill(m_frequencies, m_frequencies + 256, 0);
  std::fill(m_isVariable, m_isVariable + 256, false);
}

void Grammar::addRule(byte variable, byte first, byte second) {
  uint32 s = m_rightHandSides.size();
  m_rules.push_back(Rule(variable, s, s+2));
  m_rightHandSides.push_back(first);
  m_rightHandSides.push_back(second);
  ++m_frequencies[first];
  ++m_frequencies[second];
  ++m_frequencies[variable];
  ++m_newRules;
}

void Grammar::expandAlphabet(const std::vector<byte>& freedSymbols,
                             const std::vector<byte>& newSpecials,
                             std::vector<uint16>& nextSpecialPairs) {
  bool isNewSpecial[256] = {false};
  uint16 specialReplaces[256] = {0};

  int specialPairsLeft = specialSymbolPairsLeft();
  
  uint32 numOfSpecials = m_specialSymbols.size();
  uint32 s = 0;
  for(size_t i = 0; i < freedSymbols.size();) {
    assert(numOfSpecials*numOfSpecials >= m_specialSymbols.size());
    //if(numOfSpecials*numOfSpecials == m_specialSymbols.size()) {
    if(specialPairsLeft == 0) {
      m_specialPairReplacements.push_back(std::make_pair(false, newSpecials[s]));
      m_specialSymbols.push_back(newSpecials[s]);
      m_isSpecialSymbol[newSpecials[s]] = true;
      isNewSpecial[newSpecials[s]] = true;
      specialReplaces[newSpecials[s]] = (newSpecials[s] << 8) | newSpecials[s];
      ++s;
      ++numOfSpecials;
      specialPairsLeft = specialSymbolPairsLeft();
    } else {
      // Calculate pair for freedSymbol
      uint32 nextSpecialPair = m_specialPairReplacements.size();
      uint32 offset = (numOfSpecials -1)*(numOfSpecials -1) + 1;
      uint32 index = nextSpecialPair - offset;
      
      bool toVariable = m_isVariable[freedSymbols[i]];

      //m_specialInverse[freedSymbols[i]] = nextSpecialPair;
      m_isVariable[freedSymbols[i]] = true;
      if(index < numOfSpecials - 1) {
        /* New character is the second member of special pair
         * freedSymbol[i] maps to pair
         * (m_specialSymbols[index],m_specialSymbols.back()) */
        m_specialPairReplacements.push_back(
            std::make_pair(toVariable, freedSymbols[i]));
        isNewSpecial[freedSymbols[i]] = true;
        uint16 spPair = (m_specialSymbols[index] << 8)|m_specialSymbols.back();
        nextSpecialPairs.push_back(spPair);
        specialReplaces[freedSymbols[i]] = spPair;
      } else {
        /* New character is the first member of special pair:
         * freedSymbol[i] maps to pair
         * (m_specialSymbols.back(),m_specialSymbols[index+1-numOfSpecials]) */
        index = index+1-numOfSpecials;
        m_specialPairReplacements.push_back(
            std::make_pair(toVariable, freedSymbols[i]));
        isNewSpecial[freedSymbols[i]] = true;
        uint16 spPair = (m_specialSymbols.back() << 8)|m_specialSymbols[index];
        specialReplaces[freedSymbols[i]] = spPair;
        nextSpecialPairs.push_back(spPair);
      }
      --specialPairsLeft;
      ++i;
    }
  }
  
  assert(s<=newSpecials.size());
  std::vector<byte> nRightSides;
  uint32 leftSidesToUpdate = m_rules.size() - m_newRules;
  
  for(size_t i = 0; i < m_rules.size(); ++i) {
    if(i < leftSidesToUpdate) {
      if(!m_rules[i].isLarge() && isNewSpecial[m_rules[i].variable()]) {
        --m_frequencies[m_rules[i].variable()];
        assert(m_rules[i].variable() < 256);
        uint16 nPair = specialReplaces[m_rules[i].variable()];
        m_rules[i].changeVariable(nPair);
        ++m_frequencies[nPair >> 8];
        ++m_frequencies[nPair & 0xff];

        ++m_specialSymbolsAsVariables;
      }
    }
    uint32 posInRightSides = nRightSides.size();
    for(size_t j = m_rules[i].begin(); j < m_rules[i].end(); ++j) {
      if(isNewSpecial[m_rightHandSides[j]]) {
        --m_frequencies[m_rightHandSides[j]];
        uint16 nPair = specialReplaces[m_rightHandSides[j]];
        ++m_frequencies[nPair >> 8];
        ++m_frequencies[nPair & 0xff];
        nRightSides.push_back(nPair >> 8);
        nRightSides.push_back(nPair & 0xff);
      } else {
        nRightSides.push_back(m_rightHandSides[j]);
      }
    }
    m_rules[i].setRange(posInRightSides, nRightSides.size());
  }
  std::swap(nRightSides, m_rightHandSides);
}

void Grammar::printRules() const {
  for(size_t i = 0; i < m_rules.size(); ++i) {
    uint16 var = m_rules[i].variable();
    if(m_rules[i].isLarge()) {
      std::cout << (var >> 8) << " ";
    }
    std::cout << (var & 0xff) <<  " --> ";
    for(size_t j = m_rules[i].begin(); j < m_rules[i].end(); ++j) {
      std::cout << ((int)m_rightHandSides[j]) << " ";
    }
    std::cout << std::endl;
  }
  for(size_t i = 0; i < m_specialPairReplacements.size(); ++i) {
    if(m_specialPairReplacements[i].first) continue;
    std::cout << i << " -> " <<  ((int)m_specialPairReplacements[i].second)
              << std::endl;
  }
}

void Grammar::readGrammar(InStream* in) {
  size_t bytesRead;
  uint32 rules = utils::readPackedInteger(*in, bytesRead);
  if(rules == 0) return;
}

uint32 Grammar::writeGrammar(OutStream* out) const {
  uint32 bytes = writeNumberOfRules(out);
  if(m_rules.size() == 0) return bytes;
  bytes += writeSpecialSymbols(out);
  bytes += writeLeftSides(out);
  return bytes;
}

uint32 Grammar::writeSpecialSymbols(OutStream* out) const {
  // Number of special symbols
  out->writeByte(numberOfSpecialSymbols());
  for(std::vector<byte>::const_iterator it = m_specialSymbols.begin();
      it != m_specialSymbols.end(); ++it)
    out->writeByte(*it);
  return m_specialSymbols.size()+1;
}

uint32 Grammar::writeLeftSides(OutStream *out) const {
  uint32 bytes = writeVariableFlags(out);
  
}

uint32 Grammar::writeVariableFlags(OutStream* out) const {
  uint32 bytes = m_rules.size()/8;
  bool divisible = (m_rules.size() % 8) == 0;
  if(!divisible) ++bytes;
  byte buffer = 0;
  byte bitsLeft = 8;

  for(std::vector<Rule>::const_iterator it = m_rules.begin();
      it != m_rules.end(); ++it) {

    buffer |= ((it->isLarge()?1:0) << --bitsLeft);
    if(bitsLeft == 0) {
      bitsLeft = 8;
      out->writeByte(buffer);
      buffer = 0;
    }
  }
  if(!divisible) {
    out->writeByte(buffer);
  }
  return bytes;
}

uint32 Grammar::writeNumberOfRules(OutStream* out) const {
  int bytesNeeded;
  uint64 packedLength = utils::packInteger(m_rules.size(), &bytesNeeded);
  for(int i = 0; i < bytesNeeded; ++i) {
    out->writeByte(packedLength & 0xff);
    packedLength >>= 8;
  }
  return bytesNeeded;
}

uint32 Grammar::writeGrammar(byte* dst) const {
  uint32 s = 0;
  if(m_rules.size() > 0) {
    s += writeRightSides(dst);
    s += writeLengthsOfRules(dst + s);
  
    s += writeFreedSymbols(dst+s);

    s += writeVariables(dst + s);
    s += writeVariableFlags(dst + s);

    s += writeSpecialSymbols(dst + s);
    s += writeNumberOfSpecialSymbols(dst + s);
  }
  s += writeNumberOfRules(dst + s);
  return s;
}

void Grammar::addSpecialSymbol(byte special) {
  assert(m_updatingRules);
  m_specialSymbols.push_back(special);
}

uint32 Grammar::writeNumberOfSpecialSymbols(byte* dst) const {
  *dst = m_specialSymbols.size();
  return 1;
}

uint32 Grammar::writeSpecialSymbols(byte* dst) const {
  uint32 s = 0;
  for(std::vector<byte>::const_reverse_iterator it = m_specialSymbols.rbegin();
      it != m_specialSymbols.rend(); ++it) {
    dst[s++] = *it;
  }
  return s;
}

uint32 Grammar::writeFreedSymbols(byte* dst) const {
  uint32 s = 0;
  
  uint32 numOfFreedSymbols = 0;
  int sq = 0;
  int curr = 0;
  for(size_t i = 0; i < m_specialPairReplacements.size(); ++i) {
    if(sq == i) {
      ++curr;
      sq = curr*curr;
    } else if (!m_specialPairReplacements[i].first) {
      ++numOfFreedSymbols;
    }
  }

  bool skipped = false;
  if(numOfFreedSymbols > 0) {
    byte prev = m_specialSymbols[0];
    curr = m_specialSymbols.size()-1;
    sq = curr*curr;
    int end = 1;


    // First and second are always for special symbols
    for(int i = m_specialPairReplacements.size()-1; i > 1; --i) {
      if(sq == i) {
        --curr;
        sq = curr*curr;
        if(i == m_specialPairReplacements.size() - end) ++end;
      } else if(!m_specialPairReplacements[i].first) {
        prev = m_specialPairReplacements[i].second;
        dst[s++] = prev;
        skipped = false;
      } else {
        if(i == m_specialPairReplacements.size() - end) ++end;
        else dst[s++] = prev;
        skipped = true;
      }
    }
  }
  if(skipped) dst[s++] = 's';
  else dst[s++] = 'n';

  assert(numOfFreedSymbols <= 255);
  dst[s++] = numOfFreedSymbols;
  return s;
}

uint32 Grammar::writeNumberOfRules(byte* dst) const {
  int bytesNeeded;
  uint64 packedLength = utils::packInteger(m_rules.size(), &bytesNeeded);
  for(int i = bytesNeeded-1; i >= 0; --i) {
    dst[i] = packedLength & 0xff;
    packedLength >>= 8;
  }
  return bytesNeeded;
}

uint32 Grammar::writeVariableFlags(byte* dst) const {
  uint32 bytesNeeded = m_rules.size()/8;
  bool divisible = (m_rules.size() % 8) == 0;
  if(!divisible) ++bytesNeeded;
  uint32 currByte = bytesNeeded-1;
  byte buffer = 0;
  byte bitsInBuffer = 0;

  for(std::vector<Rule>::const_iterator it = m_rules.begin();
      it != m_rules.end(); ++it) {

    buffer |= ((it->isLarge()?1:0) << bitsInBuffer);
    ++bitsInBuffer;
    if(bitsInBuffer == 8) {
      dst[currByte--] = buffer;
      bitsInBuffer = 0;
      buffer = 0;
    }
  }
  if(!divisible) {
    dst[0] = buffer;
  }
  return bytesNeeded;
}

uint32 Grammar::writeVariables(byte* dst) const {
  uint32 s = 0;
  for(std::vector<Rule>::const_reverse_iterator it = m_rules.rbegin();
      it != m_rules.rend(); ++it) {
    uint16 var = it->variable();
    if(it->isLarge()) {
      dst[s++] = (var >> 8) & 0xff;
    }
    dst[s++] = var & 0xff;
  }
  return s;
}

/**TODO: use gamma-coding or some other more efficient code */
uint32 Grammar::writeLengthsOfRules(byte* dst) const {
  uint32 s = 0;
  for(std::vector<Rule>::const_reverse_iterator it = m_rules.rbegin();
      it != m_rules.rend(); ++it) {
    int bytesNeeded;
    uint64 packedLength = utils::packInteger((uint64)it->length(), &bytesNeeded);
    for(int i = bytesNeeded-1; i >= 0; --i) {
      dst[s+i] = packedLength & 0xff;
      packedLength >>= 8;
    }
    s += bytesNeeded;
  }
  return s;
  
}

uint32 Grammar::writeRightSides(byte* dst) const {
  uint32 s = 0;
  for(std::vector<Rule>::const_reverse_iterator it = m_rules.rbegin();
      it != m_rules.rend(); ++it) {
    std::copy(m_rightHandSides.begin() + it->begin(),
              m_rightHandSides.begin() + it->end(), dst+s);
    s += it->length();
  }
  assert(m_rightHandSides.size() == s);
  return s;
}

} //namespace bwtc
