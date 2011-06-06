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

#include "globaldefs.hpp"

using bwtc::uint64;
using bwtc::byte;
using bwtc::uint32;

/* Useful functions for debugging activites etc. */
namespace utils {

/** Floor of logarithm of base two */
byte logFloor(uint32 n);

uint32 MostSignificantBit16(uint32 n);

uint32 MostSignificantBit(uint32 n);

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
uint64 PackInteger(uint64 integer, int* bytes_needed);

uint64 UnpackInteger(uint64 packed_integer);

void WritePackedInteger(uint64 packed_integer, byte *to);

unsigned PackAndWriteInteger(uint64 integer, byte *to);

unsigned ReadAndUnpackInteger(byte *from, uint64 *to);

template <typename Integer>
void PrintBitRepresentation(Integer word) {
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
