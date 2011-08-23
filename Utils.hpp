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

#include <iostream>
#include <stack>
#include <utility>
#include <vector>

#include "globaldefs.hpp"

using bwtc::uint64;
using bwtc::byte;
using bwtc::uint32;

/* Useful functions for debugging activites etc. */
namespace utils {

/** Floor of logarithm of base two */
template <typename Unsigned>
byte logFloor(Unsigned n) {
  assert(n > 0);
  byte log = 0;
  while(n > 1) {
    n >>= 1;
    ++log;
  }
  return log;
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
  if(begin == end) {
    binaryCode(list[begin], lo, hi, bitVector);
    return;
  }
  size_t h = (end - begin) / 2;
  size_t half = begin + h;
  binaryCode(list[half], lo + h, hi + half - end, bitVector);
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

template <typename Integer, typename Input>
void binaryInterpolativeDecode(std::vector<Integer>& list, Input& input,
                               size_t lo, size_t hi, size_t elements)
{
  if(elements == 0) return;
  size_t h = (elements-1)/2;
  Integer mid = binaryDecode(input, lo + h, hi - h);
  binaryInterpolativeDecode(list, input, lo, mid-1, h);
  list.push_back(mid);
  binaryInterpolativeDecode(list, input, mid+1, hi, elements - h - 1);
}


/**Decodes binary interpolative code. It is assumed that the decoded values
 * are from range [0..maxValue].
 *
 * @param list List for the result integers.
 * @param input Source for the bits. Type must have readBit()-function
 *              returning bool.
 * @param maxValue Maximum possible value for integer to be decoded.
 * @param elements Integers to be decoded.
 */
template <typename Integer, typename Input>
void binaryInterpolativeDecode(std::vector<Integer>& list, Input& input,
                               size_t maxValue, size_t elements)
{
  binaryInterpolativeDecode(list, input, 0, maxValue, elements);
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
