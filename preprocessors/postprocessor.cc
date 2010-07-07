#include <cassert>

#include <iostream>
#include <vector>

#include "../globaldefs.h"
#include "postprocessor.h"

namespace bwtc {

uint64 UncompressCommonPairs(std::vector<byte> *compressed, uint64 length) {
  std::vector<byte>& data = *compressed;
  assert(length > 2);
  std::vector<byte> result;
  result.reserve(length - 2); /* Minimum size of result */
  /* Prepare the replacement table */
  uint16 replacements[256];
  for (unsigned i = 0; i < 256; ++i) replacements[i] = static_cast<uint16>(i);

  /* initialize value of j to the first index of compressed data */
  uint64 j = 0;
  bool escaping = false;
  byte escape_symbol;
  if (data[0] == data[1]) j = 2;
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
    std::clog << ((j - 2)/3) << " pair replacements.\n";
  }
  for( ; j < length; ++j) {
    byte current = data[j];
    if (escaping && escape_symbol == current) {
      ++j;
      result.push_back(data[j]);
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

} //namespace bwtc

