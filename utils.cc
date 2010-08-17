#include <cassert>

#include "globaldefs.h"
#include "utils.h"

namespace utils {

static const uint64 kEightBit = 0x80;

unsigned PackAndWriteInteger(uint64 integer, byte *to) {
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
unsigned ReadAndUnpackInteger(byte *from, uint64 *to) {
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

uint64 PackInteger(uint64 integer, int* bytes_needed) {
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

uint64 ReadPackedInteger(byte *to) {
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
