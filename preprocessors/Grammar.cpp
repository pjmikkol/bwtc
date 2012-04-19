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

void Grammar::increaseAlphabet(const std::vector<byte>& freedSymbols,
                               const std::vector<byte>& newSpecials,
                               std::vector<uint16>& nextSpecialPairs) {
  uint32 numOfSpecials = m_specialSymbols.size();
  uint32 s = 0;
  for(size_t i = 0; i < freedSymbols.size();) {
    if(numOfSpecials*numOfSpecials == m_specialSymbols.size()) {
      m_specialPairReplacements.push_back(std::make_pair(false, newSpecials[s]));
      m_specialSymbols.push_back(newSpecials[s]);
      m_isSpecialSymbol[newSpecials[s]] = true;
      ++s;
      ++numOfSpecials;
    } else {
      // Calculate pair for freedSymbol
      uint32 nextSpecialPair = m_specialPairReplacements.size();
      uint32 offset = (numOfSpecials -1)*(numOfSpecials -1) + 1;
      uint32 index = nextSpecialPair - offset;
      if(index < numOfSpecials -1) {
        //new character is the second member of special pair
        // freedSymbol[i] -> (m_specialSymbols[index],m_specialSymbols.back())
        m_specialPairReplacements.push_back(std::make_pair(false, freedSymbols[i]));
        nextSpecialPairs.push_back((m_specialSymbols[index] << 8) | m_specialSymbols.back());
      } else {
        //new character is the first member of special pair
        // freedSymbol[i] -> (m_specialSymbols.back(),m_specialSymbols[index+1-numOfSpecials])
        uint32 index = index+1-numOfSpecials;
        m_specialPairReplacements.push_back(std::make_pair(false, freedSymbols[i]));
        nextSpecialPairs.push_back((m_specialSymbols.back() << 8) | m_specialSymbols[index]);
      }
      ++i;
    }

  }
}


uint32 Grammar::writeGrammar(byte* dst) const {
  uint32 s = 0;
  if(m_rules.size() > 0) {
    s += writeRightSides(dst);
    s += writeLengthsOfRules(dst + s);

    s += writeFreedSymbols(dst+s);
  
    s += writeVariables(dst + s);
    s += writeLargeVariableFlags(dst + s);

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
  for(int i = 0; i < m_specialPairReplacements.size(); ++i) {
    if(sq == i) {
      ++curr;
      sq = curr*curr;
    } else if (!m_specialPairReplacements[i].first) {
      ++numOfFreedSymbols;
    }
  }
  if(numOfFreedSymbols > 0) {
    curr = m_specialSymbols.size()-1;
    sq = curr*curr;

    // First and second are always for special symbols
    for(int i = m_specialPairReplacements.size()-1; i > 1; --i) {
      if(sq == i) {
        --curr;
        sq = curr*curr;
      } else if(!m_specialPairReplacements[i].first) {
        dst[s++] = m_specialPairReplacements[i].second;
      }
    }
  }
  assert(numOfFreedSymbols <= 255);
  dst[s++] = numOfFreedSymbols;
  return s;
}

uint32 Grammar::writeNumberOfRules(byte* dst) const {
  int bytesNeeded;
  uint64 packedLength = utils::packInteger(m_rules.size(), &bytesNeeded);
  for(int i = 0; i < bytesNeeded; ++i) {
    dst[i] = packedLength & 0xff;
    packedLength >>= 8;
  }
  return bytesNeeded;
}

uint32 Grammar::writeLargeVariableFlags(byte* dst) const {
  uint32 bytesNeeded = m_rules.size()/8;
  bool divisible = (m_rules.size() % 8) != 0;
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
    for(int i = 0; i < bytesNeeded; ++i) {
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
