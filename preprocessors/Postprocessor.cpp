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

Postprocessor::Postprocessor(bool verbose) : m_verbose(verbose)
{
  std::fill(m_isSpecial, m_isSpecial + 256, false);
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

void Postprocessor::postProcess(std::vector<byte> *data) {
  PROFILE("Postprocessor::postProcess");
  uint32 grammarSize = readGrammar(&(*data)[0], data->size());
  if(grammarSize == 1) {
    data->resize(data->size()-1);
    return;
  }
  std::vector<byte>& src = *data;
  std::vector<byte> uncompressed;
  uncompressed.reserve(src.size()-1);

  uncompress(&src[0], src.size() - grammarSize, uncompressed);
  std::swap(src,uncompressed);
}

uint32 Postprocessor::readGrammar(const byte* src, size_t len) {
  const size_t last = len-1;
  //Prepare the ordinary symbols
  for(size_t i = 0; i < 256; ++i) {
    m_replacements[i].push_back((byte)i);
  }
  uint32 size = 0;
  int bytes;
  uint32 nRules = readReversedPackedInteger(src+last, &bytes);
  size += bytes;
  if(m_verbose) {
    std::clog << "Read " << nRules << " rules." << std::endl;
  }
  if(nRules == 0) return size;

  uint32 nSpecialSymbols = *(src + last - size++);
  std::vector<byte> specials;
  int specialEnumeration[256] = {0};
  for(size_t i = 0; i < nSpecialSymbols; ++i) {
    byte special = *(src + last - size++);
    specialEnumeration[special] = i;
    specials.push_back(special);
    m_isSpecial[special] = true;
  }

  // Prepare lookuptable of special pairs
  std::vector<bool> usedSpecialPair;
  usedSpecialPair.resize(nSpecialSymbols*nSpecialSymbols);
  std::fill(usedSpecialPair.begin(), usedSpecialPair.end(), false);
  for(size_t i = 0; i < nSpecialSymbols; ++i) {
    usedSpecialPair[i*i] = true;
    int replValue = (1 << 16) | (specials[i] << 8) | specials[i];
    m_replacements[replValue].push_back(specials[i]);
  }


  std::vector<bool> isLargeVariable;

  bytes = nRules/8;
  if(nRules % 8) ++bytes;
  for(int i = 0; i < bytes; ++i) {
    byte flags = *(src + last - size++);
    for(int j = 7; j >= 0; --j) {
      isLargeVariable.push_back(flags & 0x1);
      flags >>= 1;
    }
  }

  
  std::vector<std::pair<bool, uint16> > leftSides;
  for(size_t i = 0; i < nRules; ++i) {
    uint16 variable = *(src + last - size++);

    if(isLargeVariable[i]) {
      uint16 first = *(src + last - size++);
      assert(m_isSpecial[variable]);
      assert(m_isSpecial[first]);
      int s1enum = specialEnumeration[first],
          s2enum = specialEnumeration[variable];
      variable = (first << 8) | variable;
      /*
      // calculate number for the pair:
      int pairenum;
      if(s1enum > s2enum) {
        pairenum = s1enum*(s1enum+1) + s2enum + 1;
      } else {
        pairenum = s2enum*s2enum + s1enum + 1;
      }
      //usedSpecialPair[pairenum] = true;
      */
    }
    leftSides.push_back(std::make_pair(isLargeVariable[i], variable));
  }

  uint32 freedSymbols = *(src + last - size++), read = 0;
  bool off = *(src + last - size++) == 's';
  int current = 0;

  while(read < freedSymbols) {
    if(usedSpecialPair[current]) {
      ++current;
    } else {
      // Pair numbered current is used for the next symbol
      byte freedSymbol = *(src + last - size);
      ++size;
      byte next = *(src + last - size);
      if(freedSymbol == next && off && current != nSpecialSymbols*nSpecialSymbols - 1) {
        ++current;
        continue;
      }
      off = true;

      //std::cout << current << " -> " << (int)freedSymbol << std::endl; 
      
      int sqr = sqrt(current);
      int offset = sqr*sqr;
      int pairVal;
      assert(current > 1);
      if(current - offset - 1 < sqr) {
        pairVal = (specials[current - offset - 1] << 8)| specials[sqr];
      } else {
        pairVal = (specials[sqr] << 8)| specials[current-offset-sqr-1];
      }

      pairVal |= (1 << 16);
      m_replacements[pairVal].push_back(freedSymbol);
      ++current;
      ++read;
    }
  }

  std::vector<uint32> lengthsOfRules;
  for(size_t i = 0; i < nRules; ++i) {
    uint32 len = readReversedPackedInteger(src+last-size, &bytes);
    lengthsOfRules.push_back(len);
    size += bytes;
  }

  for(size_t i = 0; i < nRules; ++i) {
    int value = leftSides[i].second;
    if(leftSides[i].first) value |= (1 << 16);
    assert(lengthsOfRules[i] > 0);
    size += lengthsOfRules[i]-1;

    std::vector<byte> tmpReplacement;
    uncompress(src+last-size, lengthsOfRules[i], tmpReplacement);
    ++size;
    std::swap(m_replacements[value], tmpReplacement);

    /* Printing the rules
    if(leftSides[i].first) {
      std::cout << ((value >> 8) & 0xff) << " ";
    }
    std::cout << (value & 0xff) << " --> ";
    for(size_t j = 0; j < m_replacements[value].size(); ++j) {
      std::cout << ((int)m_replacements[value][j]) << " ";      
    }
    std::cout << std::endl;
    */
    
  }
  if(m_verbose) {
    std::clog << "Found " << nSpecialSymbols << " special symbols and "
              << freedSymbols << " freed symbols." << std::endl;
  }
  return size;
}



uint32 Postprocessor::
readReversedPackedInteger(const byte* src, int* bytesRead) {
  int bytes = 0;
  uint32 result = 0;
  while(true) {
    byte b = *src--;
    result |= ((b & 0x7f) << 7*bytes);
    ++bytes;
    if((b & 0x80) == 0) break;
  }
  *bytesRead = bytes;
  return result;
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

