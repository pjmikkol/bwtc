#include <cassert>

#include "globaldefs.h"
#include "utils.h"

namespace utils {

uint64 PackInteger(uint64 integer, int* bytes_needed) {
  /* Results for the if-clause are undefined if integer value 0x80
   * doesn't have correct type (64-bit) */
  static const uint64 kEightBit = 0x80;

  uint64 result = 0; int i;
  // For optimization (if needed) two OR-operations could be merged
  for(i = 0; integer; ++i) {
    result |= ((integer & 0x7F) << i*8);
    integer >>= 7;
    assert(i < 8);
    if (integer) result |= (kEightBit << i*8);
  }
  *bytes_needed = i;
  return result;
}

uint64 UnpackInteger(uint64 packed_integer) {
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

void WritePackedInteger(uint64 packed_integer, byte *to) {
  do {
    byte to_written = static_cast<byte>(packed_integer & 0xFF);
    packed_integer >>= 8;
    *to++ = to_written;
  } while (packed_integer);
}

} //namespace utils
