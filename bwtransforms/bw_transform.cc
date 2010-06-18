/*        Part of this code was originally released at        * 
 *        http://code.google.com/p/dcs-bwt-compressor/        */
#include "../globaldefs.h"
#include "bw_transform.h"
#include "dcbwt.h"

namespace bwtc {

byte* AllocateMemory(uint64 size) {
  return new byte[size + 1];
}

BWTransform* GiveTransformer() {
  /* When there are multiple ways to do transform this is to place to add them*/
  return new DCBWTransform();
}

/***********************************************************************
 * BwtFromSuffixArray computes the BWT from the suffix array and       *
 * writes it to result. 
 ***********************************************************************/
uint32 BwtFromSuffixArray(const byte* block, uint32 block_size,
                          const uint32* suffix_array, byte* output)
{
  byte ch = '\0';
  uint32 eob_position = 0;
  for (uint32 rank = 0; rank < block_size + 1; ++rank) {
    uint32 suffix = suffix_array[rank];
    assert(suffix <= block_size);
    if (suffix != 0) {
      ch = block[suffix - 1];
    } else {
      /* This is the end-of-block (eob) character, which cannot be
       * directly written since it has no code.
       * Instead, write a copy of the previous character here and
       * store this position in order to indicate its location. */
      eob_position = rank;
    }
    *output++ = ch;
  }
  return eob_position;
}

}
