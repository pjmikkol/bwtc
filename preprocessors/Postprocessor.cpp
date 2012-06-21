/**
 * @file Postprocessor.cpp
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
      /*std::cout << (freedSyms[i].first >> 8) << " "
                << (freedSyms[i].first & 0xff) << " -> "
                <<  ((int) freedSyms[i].second) << std::endl;*/
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


size_t Postprocessor::
uncompressCommonPairs(std::vector<byte> *compressed, size_t length) {
  static const unsigned no_repl = 70000;
  std::vector<byte>& data = *compressed;
  assert(length > 2);
  std::vector<byte> result;
  result.reserve(length - 1); /* Minimum size of result */
  /* Prepare the replacement table */
  unsigned replacements[256];
  std::fill(replacements, replacements + 256, no_repl);

  /* initialize value of j to the first index of compressed data */
  uint64 j = 0;
  bool escaping = false;
  byte escape_symbol = 0;
  int repl = (int)data[0];
  if (repl == 0) j = 1;
  else {
    unsigned i = 0;
    do {
      replacements[data[i]] = (data[i+1] << 8) | data[i+2];
      i += 3;
    } while( data[i] != data[i-3]);
    if (data[i] != data[i+1]) {
      escaping = true;
      escape_symbol = data[i+1];
    }
    j = i + 2;
  }
  if (verbosity) {
    std::clog << (!repl?0:((j - 2)/3)) << " pair replacements.";
    if(escaping) std::clog << " Escaping in use.";
    std::clog << std::endl;
  }
  for( ; j < length; ++j) {
    byte current = data[j];
    if (escaping && escape_symbol == current) {
      result.push_back(data[++j]);
    } else if (replacements[current] == no_repl) {
      result.push_back(current);
    } else {
      uint16 pair = replacements[current]&0xFFFF;
      result.push_back(pair >> 8);
      result.push_back(pair & 0xFF);
    }
  }
  data.resize(result.size());
  std::copy(result.begin(), result.end(), data.begin());
  return result.size();
}

size_t Postprocessor::
uncompressLongRuns(std::vector<byte> *compressed, size_t length) {
  std::vector<byte>& data = *compressed;
  assert(length > 2);
  std::vector<byte> result;
  result.reserve(length - 2);
  /* Runs are represented with pairs in format (length, character) */
  std::pair<unsigned, byte> replacements[256];
  for(unsigned i = 0; i < 256; ++i) 
    replacements[i] = std::make_pair(1, static_cast<byte>(i));
  bool escaping = false;
  byte escape_byte = 0;
  uint32 j; /* Init j to the start of compressed data and read replacements */
  if( data[1] == 0 ) j = 2;
  else {
    j = 0;
    byte check_val = ~data[j];
    while(1) {
      byte lengths = data[j+1];
      if(lengths == 0) {
        escape_byte = data[j];
        escaping = check_val != escape_byte;
        j += 2;
        break;
      }
      replacements[data[j]] = std::make_pair(1 << (lengths >> 4), data[j+2]);
      if ((lengths & 0x0F) == 0) {
        check_val = data[j];
        j += 4;
        escape_byte = data[j-1];
        escaping = check_val != escape_byte;
        break;
      } else {
        replacements[data[j+3]] = std::make_pair(1 << (lengths & 0x0F),
                                                 data[j+4]);
        check_val = data[j+3];
        j += 5;
      }
    }
  }
  if (verbosity) {
    unsigned repls;
    /* Compute the number of replacements from j */
    if ( (2*j-3) % 5 == 0 ) repls = (2*j-3) / 5;
    else repls = (2*j - 4) / 5;
    std::clog << repls << " run replacements. ";
    if(escaping) std::clog << "Escaping in use\n";
    else std::clog << "\n";
  }
  for(; j < length; ++j) {
    byte current = data[j];
    if(escaping && escape_byte == current) {
      result.push_back(data[++j]);
    }
    else if (replacements[current].first == 1) {
      result.push_back(current);
    }
    else {
      byte val = replacements[current].second;
      for(unsigned i = 0; i < replacements[current].first; ++i) {
        result.push_back(val);
      }
    }
  }
  data.resize(result.size());
  std::copy(result.begin(), result.end(), data.begin());
  return result.size();
}


size_t Postprocessor::
uncompressPairsAndRuns(std::vector<byte> *compressed, size_t length) {
  std::vector<byte>& data = *compressed;
  assert(length > 2);
  std::vector<byte> result;
  result.reserve(length - 3); /* Minimum size of result */
  /* Prepare the replacement table */
  Replacement replacements[256];
  std::fill(replacements, replacements + 256, Replacement(0,0,true));

  /* initialize value of j to the first index of compressed data */
  size_t j = 0;
  bool escaping = false;
  byte escapeSymbol = 0;
  if (data[0] == data[1] && data[1] == data[2]) j = 3;
  else {
    unsigned i = 0;
    do {
      replacements[data[i]] = Replacement(2, (data[i+1] << 8)|data[i+2], true);
      i += 3;
    } while( data[i] != data[i-3]);
    if (data[i] != data[i+1]) {
      escaping = true;
      escapeSymbol = data[i+1];
    }
    j = i + 2;
  }
  size_t afterPairHeader = j;
  size_t pairReplacements = (j - 2)/3;

  if(data[j+1] == 0) j += 2;
  else {
    byte checkVal = ~data[j];
    while(true) {
      byte lengths = data[j+1];
      if(lengths == 0) {
        escapeSymbol = data[j];
        escaping = checkVal != escapeSymbol;
        j += 2;
        break;
      }
      replacements[data[j]] = Replacement(1 << (lengths >> 4), data[j+2], false);
      if ((lengths & 0x0F) == 0) {
        checkVal = data[j];
        j += 4;
        escapeSymbol = data[j-1];
        escaping = checkVal != escapeSymbol;
        break;
      } else {
        replacements[data[j+3]] = Replacement(1 << (lengths & 0x0F),data[j+4],false);
        checkVal = data[j+3];
        j += 5;
      }
    }
  }
  size_t runReplacements = j - afterPairHeader;
  if((2*runReplacements-3) % 5 == 0) runReplacements = (2*runReplacements-3) /5;
  else runReplacements = (2*runReplacements-4) /5;

  if (verbosity) {
    std::clog << pairReplacements << " pair replacements and "
              << runReplacements << " run replacements.";
    if(escaping) std::clog << " Escaping in use.";
    std::clog << std::endl;
  }

  for(; j < length; ++j) {
    byte current = data[j];
    if(escaping && escapeSymbol == current) {
      result.push_back(data[++j]);
    } else if (replacements[current].length == 0) {
      result.push_back(current);
    } else if(replacements[current].isPair) {
      result.push_back(replacements[current].replacement >> 8);
      result.push_back(replacements[current].replacement & 0xff);
    } else {
      for(size_t i = 0; i < replacements[current].length; ++i)
        result.push_back(replacements[current].replacement);
    }
  }
  data.resize(result.size());
  std::copy(result.begin(), result.end(), data.begin());
  return result.size();
}

size_t Postprocessor::
uncompressSequences(std::vector<byte> *compressed, size_t length) {
  std::vector<byte>& data = *compressed;
  std::vector<byte> result;
  result.reserve(length-2);
  /* We use pair <position, length> for representing the sequences */
  std::pair<uint32, uint32> repls[256];
  std::fill(repls, repls + 256, std::make_pair(0,0));

  bool escaping = false;
  byte escapeByte = 0;
  uint64 sourcePos = 0;
  uint32 reps = 0;
  if(data[1] == 0) {
    sourcePos = 2;
  } else {
    byte prev = ~data[0];
    while(true) {
      if(prev == data[sourcePos]) break;
      prev = data[sourcePos++];
      uint64 len;
      sourcePos += utils::readAndUnpackInteger(&data[sourcePos], &len);
      repls[prev].second = len;
      repls[prev].first = sourcePos;
      sourcePos += len;
      ++reps;
    }
    escapeByte = data[++sourcePos];
    if(escapeByte != prev) escaping = true;
    ++sourcePos;
  }
  
  if(verbosity) {
    std::clog << reps << " sequences replaced. ";
    std::clog << (escaping?"U":"Not u") << "sing escape byte." << std::endl;
  }

  for(; sourcePos < length; ++sourcePos) {
    byte cur = data[sourcePos];
    uint32 len = repls[cur].second;
    if(len > 0) {
      const byte *adr = &data[repls[cur].first];
      std::copy(adr, adr + len, std::back_inserter(result));
    } else if (escaping && cur == escapeByte) {
      result.push_back(data[++sourcePos]);
    } else {
      result.push_back(data[sourcePos]);
    }
  }
  data.resize(result.size());
  std::copy(result.begin(), result.end(), data.begin());
  return data.size();
}

} //namespace bwtc

