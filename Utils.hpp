/**
 * @file Utils.hpp
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
 * Header for utility-functions which aren't related to any class that much.
 */

#ifndef BWTC_UTILS_HPP_
#define BWTC_UTILS_HPP_

#include "globaldefs.hpp"

#include <cassert>
#include <iostream>
#include <stack>
#include <utility>
#include <vector>

using bwtc::uint64;
using bwtc::byte;
using bwtc::uint32;

/* Useful functions for debugging activites etc. */
namespace utils {

static byte logFloor(unsigned n) {
  assert(n > 0);
#ifdef __GNUC__
  return  sizeof(n)*__CHAR_BIT__ - __builtin_clz(n) - 1;
#else
  byte log = 0;
  while(n > 1) {
    n >>= 1;
    ++log;
  }
  return log;
#endif
}

static byte logFloor(unsigned long n) {
  assert(n > 0);
#ifdef __GNUC__
  return  sizeof(n)*__CHAR_BIT__ - __builtin_clzl(n) - 1;
#else
  byte log = 0;
  while(n > 1) {
    n >>= 1;
    ++log;
  }
  return log;
#endif
}

/** Ceiling of logarithm of base two */
template <typename Unsigned>
byte logCeiling(Unsigned n) {
  byte log = logFloor(n);
  return static_cast<Unsigned>(1 << log) < n ? log+1 : log;
}

uint32 mostSignificantBit16(uint32 n);

uint32 mostSignificantBit(uint32 n);

/*************************************************************************
 * PackInteger and UnpackInteger                                         *
 * Packs integer to bytefields, so that the first bit of the byte        *
 * tells if the number has succeeding bytes. Few examples:               *
 * 0xF0   -> 0x01F0                                                      *
 *   -- the last byte is F0, because of the continuation-bit             *
 * 0x2    -> 0x2                                                         *
     -- no overhead here                                                 *
 * 0x142A -> 0x28AA                                                      *
 *   -- no overhead here because the most significant byte had leading   *
 *      zeroes                                                           *
 *************************************************************************/
uint64 packInteger(uint64 integer, int* bytes_needed);

uint64 unpackInteger(uint64 packed_integer);

void writePackedInteger(uint64 packed_integer, byte *to);

unsigned packAndWriteInteger(uint64 integer, byte *to);

unsigned readAndUnpackInteger(byte *from, uint64 *to);

void calculateRunFrequencies(uint64 *runFreqs, const byte *src, size_t length);

uint64 calculateRunFrequenciesAndStoreRuns(uint64 *runFreqs, byte *runseq,
  uint32 *runlen,  const byte *src, size_t length);
  
bool isPrefix(const std::string &a, const std::string &b);
  
void computeHuffmanCodes(uint32 *clen, uint32 *code);

/**Calculates the code lengths in Huffman coding. Implementation is an
 * algorithm presented in a paper "In-Place Calculation of Minimum-Redundancy
 * codes" by Alistair Moffat and Jyrki Katajainen.
 *
 * @param codeLengths Answer is returned in vector consisting of
 *                    <codelength, symbol> pairs.
 * @param freqs Array of size 256 (frequency for each character.).
 *              Array IS modified during the calculation.
 */
void calculateHuffmanLengths(std::vector<std::pair<uint64, byte> >& codeLengths,
                             uint64 *freqs);
 
template <typename Unsigned, typename BitVector>
void pushBits(Unsigned n, byte bits, BitVector& bitVector) {
  for(size_t i = 1; i <= bits; ++i) {
    bitVector.push_back((n >> (bits-i))&1);
  }
}

/**Represents given integer in bits when the lower and upper bounds for integer
 * are known.
 *
 * @param n Integer to be coded.
 * @param lo The lower bound for n (ie. n >= lo).
 * @param hi The upper bound for n (ie. n <= hi).
 * @param bits Code bits are appended into this bitvector.
 */
template <typename BitVector>
void binaryCode(size_t n, size_t lo, size_t hi, BitVector& bits) {
  size_t rangeLen = hi - lo + 1;
  if(rangeLen == 1) return;
  byte codeLength = logCeiling(rangeLen);
  size_t shortCodewords = (1 << codeLength) - rangeLen;
  size_t longCodewords2 = (rangeLen - shortCodewords)/2;
  if(n - lo < longCodewords2) {
    pushBits(n - lo, codeLength, bits);
  } else if(n - lo < longCodewords2 + shortCodewords) {
    pushBits(n - lo, codeLength - 1, bits);
  } else {
    pushBits(n - lo - shortCodewords, codeLength, bits);
  }
}

/**Calculates binary interpolative code for the sorted list of integers.
 *
 * @param list Sorted list of integers to code.
 * @param begin The first index of list to be coded.
 * @param end The last index of list to be coded.
 * @param lo Minimum possible value.
 * @param hi Maximum possible value.
 * @param bitVector The code is returned in this vector.
 */
template<typename Integer, typename BitVector>
void binaryInterpolativeCode(const std::vector<Integer>& list, size_t begin,
                             size_t end, size_t lo, size_t hi,
                             BitVector& bitVector)
{
  if(begin > end) return;
  if(end - begin == hi - lo) return;
  if(begin == end) {
    binaryCode(list[begin], lo, hi, bitVector);
    return;
  }
  size_t h = (end - begin) / 2;
  size_t half = begin + h;
  binaryCode(list[half], lo + h, hi + half - end, bitVector);
  if(half > begin)
    binaryInterpolativeCode(list, begin, half-1, lo, list[half] - 1, bitVector);
  binaryInterpolativeCode(list, half+1, end, list[half] + 1, hi, bitVector);
}

/**Calculates binary interpolative code for the sorted list of integers
 * according the paper "Binary Interpolative Coding for Effective Index
 * Compression" by Alistair Moffat and Lang Stuiver. It is assumed that
 * the coded values are from range [0..maxValue].
 *
 * @param list Sorted list of integers to code.
 * @param maxValue Maximum possible value for integer to be coded.
 * @param bitVector The code is returned in this vector.
 */
template<typename Integer, typename BitVector>
void binaryInterpolativeCode(const std::vector<Integer>& list, size_t maxValue,
                             BitVector& bitVector)
{
  assert(!list.empty());
  binaryInterpolativeCode(list, 0, list.size() - 1, 0, maxValue, bitVector);
}

template <typename Input>
size_t binaryDecode(Input& input, size_t lo, size_t hi) {
  size_t rangeLen = hi - lo + 1;
  if(rangeLen == 1) return lo;
  byte codeLength = logCeiling(rangeLen);
  size_t shortCodewords = (1 << codeLength) - rangeLen;
  size_t longCodewords2 = (rangeLen - shortCodewords)/2;
  size_t result = 0;
  for(int i = 0; i < codeLength - 1; ++i) {
    result <<= 1;
    result |= (input.readBit()) ? 1 : 0;
  }
  if(result >= longCodewords2) {
    return result + lo;
  } 
  result <<= 1;
  result |= (input.readBit()) ? 1 : 0;
  if (result < longCodewords2) return result + lo;
  else return result + lo + shortCodewords;
}

template <typename Input>
size_t binaryDecode(Input& input, size_t lo, size_t hi, size_t *bitsRead) {
  size_t rangeLen = hi - lo + 1;
  if(rangeLen == 1) return lo;
  byte codeLength = logCeiling(rangeLen);
  size_t shortCodewords = (1 << codeLength) - rangeLen;
  size_t longCodewords2 = (rangeLen - shortCodewords)/2;
  size_t result = 0;
  *bitsRead += codeLength;
  for(int i = 0; i < codeLength - 1; ++i) {
    result <<= 1;
    result |= (input.readBit()) ? 1 : 0;
  }
  if(result >= longCodewords2) {
    --*bitsRead;
    return result + lo;
  } 
  result <<= 1;
  result |= (input.readBit()) ? 1 : 0;
  if (result < longCodewords2) return result + lo;
  else return result + lo + shortCodewords;
}

template <typename Integer, typename Input>
size_t binaryInterpolativeDecode(std::vector<Integer>& list, Input& input,
                                 size_t lo, size_t hi, size_t elements)
{
  if(elements == 0) return 0;
  if(elements == hi - lo + 1) {
    size_t i = lo;
    for(Integer b = lo; i <= hi; ++i, ++b) list.push_back(b);
    return 0;
  }
  size_t bitsRead = 0;
  size_t h = (elements-1)/2;
  size_t r = elements/2 - h;
  Integer mid = binaryDecode(input, lo + h, hi - h - r, &bitsRead);

  bitsRead += binaryInterpolativeDecode(list, input, lo, mid-1, h);
  list.push_back(mid);
  bitsRead +=  binaryInterpolativeDecode(list, input, mid+1, hi, elements-h-1);
  return bitsRead;
}


/**Decodes binary interpolative code. It is assumed that the decoded values
 * are from range [0..maxValue].
 *
 * @param list List for the result integers.
 * @param input Source for the bits. Type must have readBit()-function
 *              returning bool.
 * @param maxValue Maximum possible value for integer to be decoded.
 * @param elements Integers to be decoded.
 * @return Bits read.
 */
template <typename Integer, typename Input>
size_t binaryInterpolativeDecode(std::vector<Integer>& list, Input& input,
                                 size_t maxValue, size_t elements)
{
  return binaryInterpolativeDecode(list, input, 0, maxValue, elements);
}

template <typename BitVector, typename Integer>
inline void pushBits(BitVector& bv, Integer n, size_t bits) {
  for(size_t i = 1; i <= bits; ++i)
    bv.push_back( (n >> (bits-i)) & 1 );
}

/**Forms unary code for given integer and pushes it into the back of
 * given bitvector. Few examples of unary code:
 * 1 -> 1
 * 2 -> 01
 * 7 -> 0000001
 *
 * @param to Bitvector where the result is pushed.
 * @param n Number to encode.
 */
template <typename BitVector>
inline void unaryCode(BitVector& to, size_t n) {
  while(n-- > 1) to.push_back(false);
  to.push_back(true);
}

template <typename Input>
inline size_t unaryDecode(Input& in) {
  size_t n = 1;
  while(!in.readBit()) ++n;
  return n;
}

template <typename Integer>
void printBitRepresentation(Integer word) {
  std::stack<Integer> bits;
  int bytes = sizeof(Integer);
  for(int j = 0; j < bytes; ++j) {
    for(int i = 0; i < 8; ++i) {
      int num = (word & 0x01) ? 1 : 0;
      bits.push(num);
      word >>= 1;
    }
  }
  int i = 0;
  while(!bits.empty()) {
    if(i % 8 == 0 && i) std::cout << " ";
    std::cout << bits.top();
    bits.pop();
    ++i;
  }
  std::cout << "\n";
}  
} //namespace utils

#endif
