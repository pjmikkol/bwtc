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

#include "../globaldefs.hpp"
#include "../Utils.hpp"
#include "Grammar.hpp"

#include <algorithm>
#include <cmath>

namespace bwtc {

Grammar::Grammar() : m_specialSymbolsAsVariables(0), m_updatingRules(false) {
  std::fill(m_isSpecialSymbol, m_isSpecialSymbol + 256, false);
  std::fill(m_frequencies, m_frequencies + 256, 0);
  std::fill(m_isVariable, m_isVariable + 256, false);
}

void Grammar::addRule(byte variable, byte first, byte second) {
  uint32 s = m_rightHandSides.size();
  m_rules.push_back(PrRule(variable, s, s+2));
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

/**Returns the number of special pair (first, second). */
uint32 Grammar::numberOfSpecialPair(uint32 first, uint32 second) const {
  if(first == second) {
    return first*second;
  } else if(first > second) {
    return first*(first+1) + second + 1;
  } else {
    return second*second + first + 1;
  }
}

/**Returns the pair which has ord as its number in enumeration. */
uint16 Grammar::specialPair(uint32 ord) const {
  int sqr = sqrt(ord);
  int base = sqr*sqr;
  if((int)ord == base)
    return (m_specialSymbols[sqr] << 8) | m_specialSymbols[sqr];
  int offset = ord - base - 1;
  if(offset < sqr) {
    return (m_specialSymbols[offset] << 8) | m_specialSymbols[sqr];
  } else {
    return (m_specialSymbols[sqr] << 8) | m_specialSymbols[offset-sqr];
  }
}

void Grammar::
freedSymbols(std::vector<std::pair<uint16, byte> >& replacements) const {
  assert(replacements.size() == 0);
  for(uint32 i = 0; i < m_specialPairReplacements.size(); ++i) {
    if(!m_specialPairReplacements[i].first) {
      replacements.push_back(
          std::make_pair(specialPair(i),
                         m_specialPairReplacements[i].second));
    }
  }
}

void Grammar::readGrammar(InStream* in) {
  size_t bytesRead;
  uint32 rules = utils::readPackedInteger(*in, bytesRead);
  if(rules == 0) return;
  int specialEnumeration[256] = {0};
  size_t specials = in->readByte();
  if(verbosity > 2) {
    std::clog << "Read " << specials << " special symbols." << std::endl;
  }
  // Tells if the character c (placed at index c) belongs to the original
  // alphabet (i.e. if character is part of the original string and is
  // replaced with special pair). Special characters are exception to this
  // since they can be interpreted easily enough.
  std::vector<bool> original(specials*specials, true);
  for(size_t i = 0; i < specials; ++i) {
    byte special = in->readByte();
    m_specialSymbols.push_back(special);
    m_isSpecialSymbol[special] = true;
    specialEnumeration[special] = i;
    original[i*i] = false;
  }
  uint32 maxSymbol = 0;
  
  std::vector<bool> isLargeVariable(rules);
  for(uint32 i = 0; i < rules; ++i) {
    isLargeVariable[i] = in->readBit();
  }
  in->flushBuffer();

  assert(m_rules.size() == 0);
  for(uint32 i = 0; i < rules; ++i) {
    uint16 var = in->readByte();
    ++m_frequencies[var];
    if(isLargeVariable[i]) {
      byte snd = in->readByte();
      ++m_frequencies[snd];
      uint32 specEnum = numberOfSpecialPair(
          specialEnumeration[var],specialEnumeration[snd]);
      original[specEnum] = false;
      if(specEnum > maxSymbol) maxSymbol = specEnum;
      var = (var << 8) | snd;
    } else {
      m_isVariable[var] = true;
    }
    m_rules.push_back(PrRule(var, 0, 0, isLargeVariable[i]));
  }
  
  assert(m_specialPairReplacements.size() == 0);
  size_t freedSymbols = in->readByte();
  size_t curr = 0;
  size_t sq = curr*curr;
  size_t sqr = curr;

  if(freedSymbols > 0) {
    size_t read = 0;

    while(read < freedSymbols) {
      //assert(curr < specials*specials);
      
      if(sq == curr) {
        m_specialPairReplacements.push_back(
            std::make_pair(false, m_specialSymbols[sqr]));
        ++sqr;
        sq = sqr*sqr;
      } else if(!original[curr]) {
        // special pair enumerated with curr is used as grammar variable
        // so we are not interested in its 'original' value
        m_specialPairReplacements.push_back(std::make_pair(true, 0));
      } else {
        byte orig = in->readByte();
        m_specialPairReplacements.push_back(std::make_pair(false, orig));
        ++read;
      }
      ++curr;
    }
  }
  for(;curr <= maxSymbol; ++curr) {
    if(sq == curr && sqr < m_specialSymbols.size()) {
      m_specialPairReplacements.push_back(
          std::make_pair(false, m_specialSymbols[sqr]));
      ++sqr;
      sq = sqr*sqr;
    } else {
      m_specialPairReplacements.push_back(std::make_pair(true, 0));
    }
  }

  //Assume that the lengths are between 2 and 4
  size_t rightSidesLength = 0;
  size_t add = (rules%4 == 0)?0:1;
  for(size_t i = 0; i < rules/4 + add; ++i) {
    byte lb = in->readByte();
    for(size_t j = 0; j < 4; ++j) {
      size_t ruleIndex = 4*i + j;
      if(ruleIndex >= m_rules.size()) continue;
      size_t l = 2 + ((lb >> (6 - 2*j)) & 0x3);
      m_rules[ruleIndex].setRange(rightSidesLength, rightSidesLength+l);
      rightSidesLength += l;
    }
  }
  m_rightHandSides.resize(rightSidesLength);

  for(size_t i = 0; i < rightSidesLength; ++i) {
    m_rightHandSides[i] = in->readByte();
    ++m_frequencies[m_rightHandSides[i]];
  }
}

uint32 Grammar::writeGrammar(OutStream* out) const {
  uint32 bytes = writeNumberOfRules(out);
  if(m_rules.size() == 0) return bytes;
  bytes += writeSpecialSymbols(out);
  bytes += writeLeftSides(out);
  bytes += writeFreedSymbols(out);
  bytes += writeLengthsOfRules(out);
  bytes += writeRightSides(out);

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
  return bytes + writeVariables(out);
}

uint32 Grammar::writeVariables(OutStream* out) const {
  uint32 bytes = 0;
  for(std::vector<PrRule>::const_iterator it = m_rules.begin();
      it != m_rules.end(); ++it) {
    uint16 var = it->variable();
    if(it->isLarge()) {
      out->writeByte((var >> 8) & 0xff);
      ++bytes;
    }
    out->writeByte(var & 0xff);
    ++bytes;
  }
  return bytes;
}

uint32 Grammar::numberOfFreedSymbols() const {
  uint32 numOfFreedSymbols = 0;
  size_t sq = 0;
  size_t curr = 0;
  for(size_t i = 0; i < m_specialPairReplacements.size(); ++i) {
    if(sq == i) {
      ++curr;
      sq = curr*curr;
    } else if (!m_specialPairReplacements[i].first) {
      ++numOfFreedSymbols;
    }
  }
  return numOfFreedSymbols;
}

uint32 Grammar::writeFreedSymbols(OutStream *out) const {
  uint32 freedSymbols = numberOfFreedSymbols();
  assert(freedSymbols <= 255);
  out->writeByte(freedSymbols);
  uint32 bytes = 1;
  
  if(freedSymbols > 0) {
    int curr = 2;
    int sq = curr*curr;
    int sqr = curr;

    while(freedSymbols > 0) {
      if(sq == curr) {
        ++sqr;
        sq = sqr*sqr;
      } else if(!m_specialPairReplacements[curr].first) {
        byte symbol = m_specialPairReplacements[curr].second;
        out->writeByte(symbol);
        --freedSymbols;
        ++bytes;
      }
      ++curr;
    }
  }
  return bytes;
}

/**Assume that rules have length of 2 to 4 (PairReplacer satisfies this).*/
uint32 Grammar::writeLengthsOfRules(OutStream* out) const {
  byte buffer = 0;
  uint32 rulesInBuffer = 0;
  for(std::vector<PrRule>::const_iterator it = m_rules.begin();
      it != m_rules.end(); ++it) {
    assert(it->length() <= 4 && it->length() >= 2);
    buffer <<= 2;
    //0 -> 2, 1 -> 3, 2 -> 4
    buffer |= (it->length() - 2);
    if(++rulesInBuffer == 4) {
      out->writeByte(buffer);
      rulesInBuffer = 0;
      buffer = 0;
    }
  }
  uint32 bytes = m_rules.size()/4;
  if(rulesInBuffer != 0) {
    ++bytes;
    out->writeByte(buffer << (8 - 2*rulesInBuffer));
  }
  return bytes;
}

uint32 Grammar::writeRightSides(OutStream* out) const {
  for(std::vector<byte>::const_iterator it = m_rightHandSides.begin();
      it != m_rightHandSides.end(); ++it) {
    out->writeByte(*it);
  }
  return m_rightHandSides.size();
}

uint32 Grammar::writeVariableFlags(OutStream* out) const {
  uint32 bytes = m_rules.size()/8;
  bool divisible = (m_rules.size() % 8) == 0;
  if(!divisible) ++bytes;
  byte buffer = 0;
  byte bitsLeft = 8;

  for(std::vector<PrRule>::const_iterator it = m_rules.begin();
      it != m_rules.end(); ++it) {

    buffer |= ((it->isLarge()?1:0) << --bitsLeft);
    if(bitsLeft == 0) {
      out->writeByte(buffer);
      bitsLeft = 8;
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

void Grammar::addSpecialSymbol(byte special) {
  assert(m_updatingRules);
  m_specialSymbols.push_back(special);
}

} //namespace bwtc
