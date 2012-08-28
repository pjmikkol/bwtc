/**
 * @file Postprocessor.cpp
 * @author Pekka Mikkola <pmikkol@gmail.com>
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
 * Implementations for the inverses of preprocessing algorithms.
 */

#include <cassert>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iterator>
#include <map>
#include <utility>
#include <vector>

#include "../globaldefs.hpp"
#include "Postprocessor.hpp"
#include "../Utils.hpp"
#include "../Profiling.hpp"
#include "Grammar.hpp"

namespace bwtc {

Postprocessor::Replacement::Replacement()
    : length(0), replacement(0), isPair(false) {}

Postprocessor::Replacement::Replacement(uint32 len, uint16 repl, bool pair)
    : length(len), replacement(repl), isPair(pair) {}

Postprocessor::Replacement::Replacement(const Replacement& r)
    : length(r.length), replacement(r.replacement), isPair(r.isPair) {}

Postprocessor::Replacement&
Postprocessor::Replacement::operator=(const Postprocessor::Replacement& r) {
  length = r.length;
  replacement = r.replacement;
  isPair = r.isPair;
  return *this;
}

Postprocessor::Postprocessor(bool verbose, const Grammar& grammar)
    : m_verbose(verbose), m_hasRules(false) {
  for(size_t i = 0; i < 256; ++i) {
    m_replacements[i].push_back((byte)i);
  }

  for(size_t i = 0; i < 256; ++i) {
    m_isSpecial[i] = grammar.isSpecial(i);
  }
  if(grammar.numOfRules() > 0) m_hasRules = true;
  
  // Freed symbols
  const int highBit = 1 << 16;
  {
    std::vector<std::pair<uint16, byte> > freedSyms;
    grammar.freedSymbols(freedSyms);
    for(size_t i = 0; i < freedSyms.size(); ++i) {
      assert(m_replacements[highBit | freedSyms[i].first].size() == 0);
      m_replacements[highBit | freedSyms[i].first].
          push_back(freedSyms[i].second);
    }
  }

  // Rules
  for(size_t i = 0; i < grammar.numberOfRules(); ++i) {
    Grammar::Rule rule = grammar.getRule(i);
    std::vector<byte> tmp;
    uncompress(rule.begin(), rule.length(), tmp);
    int repAddress = rule.variable();
    if(rule.isLarge())  repAddress |= highBit;
    std::swap(m_replacements[repAddress], tmp);
  }
}



void Postprocessor::
uncompress(const byte* src, size_t length, std::vector<byte>& dst) {
  for(size_t i = 0; i < length; ++i) {
    int key = src[i];
    if(m_isSpecial[src[i]]) {
      key = (1 << 16) | (src[i] << 8 )| src[i+1];
      ++i;
    }
    for(size_t j = 0; j < m_replacements[key].size(); ++j) {
      dst.push_back(m_replacements[key][j]);
    }
  }
}

size_t Postprocessor::
uncompress(const byte* data, size_t length, OutStream* to) const {
  PROFILE("Postprocessor::uncompress");
  if(!m_hasRules) {
    to->writeBlock(data, data+length);
    return length;
  }
  
  size_t uncompressedLength = 0;
  for(size_t i = 0; i < length; ++i) {
    int key = data[i];
    if(m_isSpecial[key]) {
      key = (1 << 16) | (data[i] << 8 )| data[i+1];
      ++i;
    }
    size_t len = m_replacements[key].size();
    to->writeBlock(&m_replacements[key][0], &m_replacements[key][len]);
    uncompressedLength += len;
  }
  return uncompressedLength;
}

} //namespace bwtc

