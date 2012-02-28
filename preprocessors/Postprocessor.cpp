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
#include <iostream>
#include <iterator>
#include <map>
#include <utility>
#include <vector>

#include "../globaldefs.hpp"
#include "Postprocessor.hpp"
#include "../Utils.hpp"

namespace bwtc {

PostProcessor::Replacement::Replacement()
    : length(0), replacement(0), isPair(false) {}

PostProcessor::Replacement::Replacement(uint32 len, uint16 repl, bool pair)
    : length(len), replacement(repl), isPair(pair) {}

PostProcessor::Replacement::Replacement(const Replacement& r)
    : length(r.length), replacement(r.replacement), isPair(r.isPair) {}

PostProcessor::Replacement&
PostProcessor::Replacement::operator=(const PostProcessor::Replacement& r) {
  length = r.length;
  replacement = r.replacement;
  isPair = r.isPair;
  return *this;
}

PostProcessor::PostProcessor(const std::string& postProcOptions)
    : m_options(postProcOptions) {}

void PostProcessor::postProcess(std::vector<byte> *data) {
  // TODO: optimize unnecessary copying away
  size_t length = data->size();
  for(size_t i = 0; i < m_options.size(); ++i) {
    if(m_options[i] == 'r') {
      length = uncompressLongRuns(data, length);
    } else if (m_options[i] == 'p') {
      length = uncompressCommonPairs(data, length);
    } else if (m_options[i] == 'c') {
      length = uncompressPairsAndRuns(data, length);
    } else if (m_options[i] == 's') {
      length = uncompressSequences(data, length);
    }
  }
}


size_t PostProcessor::
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

size_t PostProcessor::
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


size_t PostProcessor::
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

size_t PostProcessor::
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

