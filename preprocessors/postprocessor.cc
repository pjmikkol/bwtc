/**************************************************************************
 *  Copyright 2010, Pekka Mikkola, pjmikkol (at) cs.helsinki.fi           *
 *                                                                        *
 *  This file is part of bwtc.                                            *
 *                                                                        *
 *  bwtc is free software: you can redistribute it and/or modify          *
 *  it under the terms of the GNU General Public License as published by  *
 *  the Free Software Foundation, either version 3 of the License, or     *
 *  (at your option) any later version.                                   *
 *                                                                        *
 *  bwtc is distributed in the hope that it will be useful,               *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *  GNU General Public License for more details.                          *
 *                                                                        *
 *  You should have received a copy of the GNU General Public License     *
 *  along with bwtc.  If not, see <http://www.gnu.org/licenses/>.         *
 **************************************************************************/

#include <cassert>

#include <iostream>
#include <map>
#include <utility>
#include <vector>

#include "../globaldefs.h"
#include "postprocessor.h"
#include "../utils.h"

namespace bwtc {

uint64 UncompressCommonPairs(std::vector<byte> *compressed, uint64 length) {
  std::vector<byte>& data = *compressed;
  assert(length > 2);
  std::vector<byte> result;
  //result.reserve(length - 3); /* Minimum size of result */
  /* Prepare the replacement table */
  uint16 replacements[256];
  for (unsigned i = 0; i < 256; ++i) replacements[i] = static_cast<uint16>(i);

  /* initialize value of j to the first index of compressed data */
  uint64 j = 0;
  bool escaping = false;
  byte escape_symbol;
  if (data[0] == data[1] && data[1] == data[2]) j = 3;
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
    std::clog << ((j - 2)/3) << " pair replacements.";
    if(escaping) std::clog << " Escaping in use.";
    std::clog << "\n";
  }
  for( ; j < length; ++j) {
    byte current = data[j];
    if (escaping && escape_symbol == current) {
      result.push_back(data[++j]);
    }
    else if (replacements[current] == current) {
      result.push_back(current);
    }
    else {
      uint16 pair = replacements[current];
      result.push_back(pair >> 8);
      result.push_back(pair & 0xFF);
    }
  }
  data.resize(result.size());
  std::copy(result.begin(), result.end(), data.begin());
  return result.size();
}

uint64 UncompressLongRuns(std::vector<byte> *compressed, uint64 length) {
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
  uint64 j; /* Init j to the start of compressed data and read replacements */
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
  if (verbosity > 2) {
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

uint64 UncompressSequences(std::vector<byte> *compressed, uint64 length) {
  std::vector<byte>& data = *compressed;
  std::vector<byte> result;
  result.reserve(length/2);
  /* We use pair <position, length> for representing the sequences */
  std::pair<uint64, unsigned> repls[256];
  for(int i = 0; i < 256; ++i) {
    repls[i] = std::pair<uint64, unsigned>(0,1);
  }
  bool escaping = false;
  byte escape_byte;
  uint64 source_pos = 0;
  if(data[1] == 0) {
    source_pos = 2;
  } else {
    byte prev = ~data[0];
    while(1) {
      if(prev == data[source_pos]) break;
      prev = data[source_pos++];
      uint64 len;
      source_pos += utils::ReadAndUnpackInteger(&data[source_pos], &len);
      repls[prev].second = len;
      repls[prev].first = source_pos;
      source_pos += len;
    }
    escape_byte = data[++source_pos];
    if(escape_byte != prev) escaping = true;
    ++source_pos;
  }

  if(verbosity > 2) {
    int rep = 0;
    for(int i = 0; i < 256; ++i) {
      if(repls[i].second > 1) ++rep;
    }
    std::clog << rep << " sequences replaced. ";
    if(escaping) std::clog << "U";
    else std::clog << "Not u";
    std::clog << "sing escape byte.\n";
  }

  for(; source_pos < length; ++source_pos) {
    if(repls[data[source_pos]].second > 1) {
      for(unsigned i = 0; i < repls[data[source_pos]].second; ++i)
        result.push_back(data[repls[data[source_pos]].first + i]);
    } else if (escaping && data[source_pos] == escape_byte) {
      result.push_back(data[++source_pos]);
    } else {
      result.push_back(data[source_pos]);
    }
  }
  data.resize(result.size());
  std::copy(result.begin(), result.end(), data.begin());
  return data.size();
}


} //namespace bwtc

