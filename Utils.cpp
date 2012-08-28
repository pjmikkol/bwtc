/**
 * @file Utils.cpp
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
 * Implementation of utility functions.
 */

#include <cassert>
#include <algorithm>
#include <utility>
#include <map>
#include <vector>

#include "globaldefs.hpp"
#include "Utils.hpp"

namespace utils {

uint32 mostSignificantBit16(uint32 n) {
  assert(n < (1 << 16));
  n |= (n >> 1);
  n |= (n >> 2);
  n |= (n >> 4);
  n |= (n >> 8);
  return n & (~(n >> 1));
}

uint32 mostSignificantBit(uint32 n) {
  assert(sizeof(uint32) == 4);
  n |= (n >> 1);
  n |= (n >> 2);
  n |= (n >> 4);
  n |= (n >> 8);
  n |= (n >> 16);
  return n & (~(n >> 1));
}

static const uint64 kEightBit = 0x80;

unsigned packAndWriteInteger(uint64 integer, byte *to) {
  unsigned bytes_used = 0;
  while(1) {
    if(integer) {
      byte to_written = static_cast<byte>(integer & 0x7F);
      integer >>= 7;
      if(integer) to_written |= 0x80;
      *to++ = to_written;
      ++bytes_used;
    } else {
      return bytes_used;
    }
  }
}

/* Returns the bytes used for the integer */
unsigned readAndUnpackInteger(byte *from, uint64 *to) {
  unsigned bytes_read = 0;
  uint64 result = 0;
  while(1) {
    uint64 read = *from;
    result |= ((read & 0x7F) << 7*bytes_read);
    ++bytes_read;
    if (((*from) & 0x80) == 0) break;
    else ++from;
  }
  *to = result;
  return bytes_read;
}

uint64 packInteger(uint64 integer, int* bytes_needed) {
  uint64 result = 0; int i;
  // For optimization (if needed) two OR-operations could be merged
  for(i = 0; integer; ++i) {
    result |= ((integer & 0x7F) << i*8);
    integer >>= 7;
    assert(i < 8);
    if (integer) result |= (kEightBit << i*8);
  }
  if(i == 0) ++i;
  *bytes_needed = i;
  return result;
}

uint64 unpackInteger(uint64 packed_integer) {
  uint64 result = 0; int bits_handled = 0;
  bool bits_left;
  do {
    bits_left = (packed_integer & 0x80) != 0;
    result |= ((packed_integer & 0x7F) << bits_handled);
    packed_integer >>= 8;
    bits_handled += 7;
    assert(bits_handled <= 56);
  } while(bits_left);
  return result;
}

void calculateRunFrequencies(uint64 *runFreqs, const byte *src, size_t length)
{
  const byte *prev = src;
  const byte *curr = src+1;
  do {
    while(curr < src + length && *prev == *curr) ++curr;
    ++runFreqs[*prev];
    prev = curr++;
  } while(curr < src + length);
  if(prev < src + length) ++runFreqs[*prev];
}

size_t calculateRunsAndCharacters(uint64 *runFreqs, const byte *src,
                                  size_t length, std::map<uint32, uint32>& runs)
{
  size_t totalRuns = 0;
  const byte *prev = src;
  const byte *curr = src+1;
  do {
    while(curr < src + length && *prev == *curr) ++curr;
    ++totalRuns;
    ++runFreqs[*prev];
    ++runs[curr - prev];
    prev = curr++;
  } while(curr < src + length);
  if(prev < src + length) {
    ++runFreqs[*prev];
    ++runs[1];
    ++totalRuns;
  }
  return totalRuns;
}


uint64 calculateRunFrequenciesAndStoreRuns(uint64 *runFreqs, byte *runseq,
  uint32 *runlen, const byte *src, size_t length)
{
  uint64 runs_cnt = 0;
  byte *ptr = runseq;
  const byte *prev = src;
  const byte *curr = src+1;
  do {
    while(curr < src + length && *prev == *curr) ++curr;
    ++runFreqs[*prev];
    *ptr++ = *prev;
    runlen[runs_cnt++] = curr - prev;
    prev = curr++;
  } while(curr < src + length);
  if(prev < src + length) {
    ++runFreqs[*prev];
    runlen[runs_cnt++] = 1;
    *ptr++ = *prev;
  }
  return runs_cnt;
}

bool isPrefix(const std::string &a, const std::string &b) {
  if (a.length() > b.length()) return false;
  bool isPref = true;
  for (int j = 0; j < (int)a.length(); ++j)
    if (a[j] != b[j]) isPref = false;
  return isPref;
}

void computeHuffmanCodes(uint32 *clen, uint32 *code) {
  // Compute the number of codes with given length.
  uint32 lengthsCount[256];
  std::fill(lengthsCount, lengthsCount + 256, 0);
  for (bwtc::int32 k = 0; k < 256; ++k)
    ++lengthsCount[clen[k]];
  bwtc::int32 maxCodeLen = 0;
  for (bwtc::int32 k = 0; k < 256; ++k)
    if (lengthsCount[k] > 0) maxCodeLen = k;

  // Compute codes. Only 8bit-suffix of each code is
  // stored. The remaining bits are quaranteed to be 0.
  uint32 start[256];
  std::fill(start, start + 256, 0);
  uint32 first = 0x00;
  for (bwtc::int32 k = maxCodeLen; k >= 0; --k) {
    start[k] = first;
    first = (first + lengthsCount[k]) >> 1;
  }
  std::fill(code, code + 256, 0);
  for (bwtc::int32 k = 0; k < 256; ++k)
    if (clen[k] > 0) code[k] = start[clen[k]]++;

  // Check if the code is prefix free.
  /*std::string full_codes[256];
  for (int32 k = 0; k < 256; ++k)
    full_codes[k] = "";
  for (int32 k = 0; k < 256; ++k) {
    if (clen[k] > 8) {
      for (int32 rr = 0; rr < (int32)clen[k] - 8; ++rr)
        full_codes[k] += "0";
    }
    int32 limit = std::min(static_cast<uint32>(8), clen[k]);
    for (int32 rr = 0; rr < limit; ++rr)
      if (code[k] & (1 << (limit - 1 - rr))) full_codes[k] += "1";
      else full_codes[k] += "0";
  }
  for (int32 k = 0; k < 256; ++k) {
    for (int32 p = 0; p < 256; ++p) {
      if (k != p && clen[k] > 0 && clen[p] > 0) {
        if (isPrefix(full_codes[k], full_codes[p]) ||
            isPrefix(full_codes[p], full_codes[k])) {
          fprintf(stderr,"ERROR.\n");
          fprintf(stderr,"full_codes[%d] = %s\n", k, full_codes[k].c_str());
          fprintf(stderr,"full_codes[%d] = %s\n", p, full_codes[p].c_str());
          fprintf(stderr,"clen:\n");
          for (int32 t = 0; t < 256; ++t)
            fprintf(stderr,"  clen[%d] = %u\n", t, clen[t]);
          exit(1);
        }
      }
    }
  }*/
}

namespace {

struct HuTuckerNode {
  HuTuckerNode(HuTuckerNode* l, HuTuckerNode* r, uint64 w, uint32 s)
      : left(l), right(r), weight(w), symbol(s) {}

  HuTuckerNode* left;
  HuTuckerNode* right;
  uint64 weight;
  uint32 symbol;
};

int nextTerminal(const std::vector<HuTuckerNode*>& v, int i) {
  int s = v.size();
  while(i < s && (!v[i] || v[i]->left)) ++i;
  if (i == s) --i;
  return i;
}

int nextValid(const std::vector<HuTuckerNode*>& v, int i) {
  int s = v.size(); 
  while(i < s && !v[i]) ++i;
  return i;
}

std::pair<int, int> findMinSum(const std::vector<HuTuckerNode*>& v, int begin,
                               int end, uint64& minSum) {
  int s = nextValid(v, begin+1);
  int left = begin, right = s;
  uint64 mS = v[s]->weight + v[begin]->weight;
  while(s <= end) {
    uint64 sum = v[s]->weight + v[begin]->weight;
    if(sum < mS) {
      mS = sum; left = begin, right = s;
    }

    if(v[s]->weight < v[begin]->weight) {
      begin = s;
    }
    s = nextValid(v, s+1);
  }

  minSum = mS;
  return std::make_pair(left, right);
}

void collectLengthsAndDestroyTree(HuTuckerNode* node, uint32 d,
                                  std::map<uint32, uint32>& lengths) {
  if(node->left == 0) {
    lengths[node->symbol] = d;
    delete node;
  } else {
    collectLengthsAndDestroyTree(node->left, d+1, lengths);
    collectLengthsAndDestroyTree(node->right, d+1, lengths);
    delete node;
  }
}
} // namespace

void calculateHuTuckerLengths(
    std::vector<std::pair<uint64, uint32> >& codeLengths,
    uint64 *freqs, const std::vector<uint32>& names)
{
  if(names.size() == 1) {
    codeLengths.push_back(std::make_pair(1,names[0]));
    return;
  } else if(names.size() == 2) {
    codeLengths.push_back(std::make_pair(1,names[0]));
    codeLengths.push_back(std::make_pair(1,names[1]));
    return;
  }

  std::vector<HuTuckerNode*> nodes;
  for(size_t i = 0; i < names.size(); ++i) {
    nodes.push_back(new HuTuckerNode(0, 0, freqs[i], names[i]));
  }
  size_t removed = 0;
  size_t im1 = 0, curr = 0;
  // first pass to combine nodes at the beginning:
  if(nodes[0]->weight <= nodes[2]->weight) {
    HuTuckerNode* n = new HuTuckerNode(
        nodes[0], nodes[1], nodes[0]->weight +  nodes[1]->weight, 0);
    nodes[0] = n;
    nodes[1] = 0;
    curr = 2;
    ++removed;
  } else {
    curr = 1;
  }
  while(curr + 2 < nodes.size()) {
    if(nodes[im1]->weight > nodes[curr+1]->weight &&
       nodes[curr]->weight <= nodes[curr+2]->weight)
    {
      HuTuckerNode* n = new HuTuckerNode(
          nodes[curr], nodes[curr + 1],
          nodes[curr]->weight + nodes[curr + 1]->weight, 0);
      nodes[curr] = n;
      nodes[curr+1] = 0;
      im1 = curr;
      curr += 2;
      ++removed;
    } else {
      im1 = curr;
      ++curr;
    }
  }
  if(curr == nodes.size() - 2) {
    if(nodes[curr + 1]->weight < nodes[curr - 1]->weight) {
      HuTuckerNode* n = new HuTuckerNode(
          nodes[curr], nodes[curr + 1],
          nodes[curr]->weight + nodes[curr + 1]->weight, 0);
      nodes[curr] = n;
      nodes[curr+1] = 0;
      ++removed;
    }
  }

  size_t deleted = removed;
  while(removed < names.size()-1) {
    if(deleted >= nodes.size()/2) {
      int i = 0;
      for(size_t j = 0; j < nodes.size(); ++j) {
        if(nodes[j]) nodes[i++] = nodes[j];
      }
      nodes.resize(i);
      deleted = 0;
    }
    int first = nextTerminal(nodes, 0);
    int second = nextTerminal(nodes, first+1);

    uint64 minWeight = nodes[first]->weight + nodes[second]->weight;
    int minLeft = first, minRight = second;

    while(first + 1 < (int)nodes.size()) {
      uint64 minSum;
      std::pair<int, int> lr = findMinSum(nodes, first, second, minSum);
      if(minSum < minWeight) {
        minWeight = minSum;
        minLeft = lr.first, minRight = lr.second;
      }
      first = second;
      second = nextTerminal(nodes, first+1);
    }

    HuTuckerNode* n = new HuTuckerNode(
        nodes[minLeft], nodes[minRight], minWeight, 0);
    nodes[minLeft] = n;
    nodes[minRight] = 0;
    
    ++removed;
    ++deleted;
  }

  std::map<uint32, uint32> lengths;
  collectLengthsAndDestroyTree(nodes[0], 0, lengths);

  for(std::map<uint32, uint32>::const_iterator it = lengths.begin();
      it != lengths.end(); ++it)
  {
    std::cout << it->first << " " << it->second << std::endl;
    codeLengths.push_back(std::make_pair(it->second, it->first));
  }
}

void calculateHuffmanLengths(std::vector<std::pair<uint64, uint32> >& codeLengths,
                             uint64 *freqs, const std::vector<uint32>& names)
{
  assert(codeLengths.size() == 0);
  for(size_t i = 0; i < names.size(); ++i) {
    if(freqs[i])
      codeLengths.push_back(std::make_pair(freqs[i], names[i]));
  }
  calculateCodeLengths(codeLengths, freqs);
}

    
void calculateHuffmanLengths(std::vector<std::pair<uint64, uint32> >& codeLengths,
                             uint64 *freqs, uint32 alphabetSize)
{
  assert(codeLengths.size() == 0);
  for(size_t i = 0; i < alphabetSize; ++i) {
    if(freqs[i])
      codeLengths.push_back(std::make_pair(freqs[i], i));
  }
  calculateCodeLengths(codeLengths, freqs);
}

void calculateCodeLengths(std::vector<std::pair<uint64, uint32> >& codeLengths,
                          uint64 *freqs, bool sorted)
{
  assert(codeLengths.size() >  0);
  if(codeLengths.size() == 1) {
    codeLengths[0].first = 1;
    return;
  }
  bool allocate = !freqs;
  if(allocate) freqs = new uint64[codeLengths.size()];

  if(!sorted) std::sort(codeLengths.begin(), codeLengths.end());
  const size_t n = codeLengths.size();
  for(size_t i = 0; i < n; ++i) {
    freqs[i] = codeLengths[i].first;
  }

  size_t s = 0, r = 0;
  for(size_t t = 0; t < n-1; ++t) {
    if(s >= n || (r < t && freqs[r] < freqs[s])) {
      freqs[t] = freqs[r];
      freqs[r++] = t;
    } else {
      freqs[t] = freqs[s++];
    }
    if(s >= n || (r < t && freqs[r] < freqs[s])) {
      freqs[t] += freqs[r];
      freqs[r++] = t;
    } else {
      freqs[t] += freqs[s++];
    }
  }

  freqs[n-2] = 0;
  for(int t = n-3; t >= 0; --t) {
    freqs[t] = freqs[freqs[t]] + 1;
  }
  int a = 1, u = 0;
  size_t depth = 0;
  int x = n-1, t = n-2;
  while(a > 0) {
    while(t >= 0 && freqs[t] == depth) { ++u; --t; }
    while(a > u) { freqs[x] = depth; --x; --a; }
    a = 2*u;
    ++depth;
    u = 0;
  }
  for(size_t i = 0; i < n; ++i) {
    codeLengths[i].first = freqs[i];    
  }
  if(allocate) delete [] freqs;
}

void writePackedInteger(uint64 packed_integer, byte *to) {
  do {
    byte to_written = static_cast<byte>(packed_integer & 0xFF);
    packed_integer >>= 8;
    *to++ = to_written;
  } while (packed_integer);
}

uint64 readPackedInteger(byte *to) {
  static const uint64 kEndSymbol = static_cast<uint64>(1) << 63;
  static const uint64 kEndMask = static_cast<uint64>(1) << 7;

  uint64 packed_integer = 0;
  bool bits_left = true;
  int i;
  for(i = 0; bits_left; ++i) {
    uint64 read = *to++;
    bits_left = (read & kEndMask) != 0;
    packed_integer |= (read << i*8);
  }
  if (packed_integer == 0x80) return kEndSymbol;
  return packed_integer;

}

} //namespace utils
