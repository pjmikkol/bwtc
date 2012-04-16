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
#include "Grammar.hpp"

namespace bwtc {

Grammar::Grammar() : m_usedSpecialSymbolPairs(0) {
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

} //namespace bwtc
