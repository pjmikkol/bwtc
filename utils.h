#ifndef BWTC_UTILS_H_
#define BWTC_UTILS_H_

#include <iostream>
#include <stack>

#include "globaldefs.h"

/* Useful functions for debugging activites etc. */
namespace utils {

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
