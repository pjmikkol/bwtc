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

Grammar::Grammar() : m_specialSymbolsAsVariables(0) {
  std::fill(m_isSpecialSymbol, m_isSpecialSymbol + 256, false);
  std::fill(m_frequencies, m_frequencies + 256, 0);
}

void Grammar::addRule(byte variable, byte first, byte second) {
  uint32 s = m_rightHandSides.size();
  m_rules.push_back(Rule(variable, s, s+2));
  m_rightHandSides.push_back(first);
  m_rightHandSides.push_back(second);
}

void Grammar::addRule(byte variable, byte* begin, size_t length) {
  uint32 s = m_rightHandSides.size();
  m_rules.push_back(Rule(variable, s, s+length));
  std::copy(begin, begin+length, std::back_inserter(m_rightHandSides));
}

void Grammar::writeGrammar(byte* dst) const {
  uint32 s = 0;
  if(m_rules.size() > 0) {
    s += writeRightSides(dst);
    s += writeLengthsOfRules(dst + s);
    s += writeVariables(dst + s);
    s += writeLargeVariableFlags(dst + s);
  }
  s += writeNumberOfRules(dst + s);
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
