#include <cassert>

#include <numeric>  // for partial_sum
#include <vector>

#include "../globaldefs.h"
#include "inverse_bwt.h"

/*************************************************************
 * Parts of this code have been previously released here:    *
 * http://code.google.com/p/dcs-bwt-compressor/              *
 *************************************************************/

namespace bwtc {

InverseBWTransform* GiveInverseTransformer() {
  return new FastInverseBWTransform();
}

std::vector<byte>* InverseBWTransform::AllocateMemory(uint64 block_size) {
  return new std::vector<byte>(block_size);
}

uint64 FastInverseBWTransform::MaxBlockSize(uint64 memory_budget) const {
  memory_budget -= kMemoryOverhead;
  return memory_budget / sizeof(uint32);
}

std::vector<byte>* FastInverseBWTransform::DoTransform(const byte* bwt,
                                                       uint64 bwt_size,
                                                       uint64 eob_position)
{
  // rank[i] will be the number of occurrences of bwt[i] in bwt[0..i-1]
  std::vector<uint32> bwt_rank_low24(bwt_size);
  std::vector<uint32> rank_milestone_buffer[256];

  // count[] serves two purposes:
  // 1. When the scan of the BWT reaches position i,
  //    count[c+1] is the number of occurrences of c in bwt[0..i-1]
  //    (count[0] counts the EOB symbol).
  // 2. During the ouput generation, count[c] is the total number
  //    of characters smaller than c (including the single EOB symbol).
  std::vector<uint32> count(257, 0);

  // count EOB
  bwt_rank_low24[eob_position] = 0;
  count[0] = 1;
  // count other characters
  for (uint32 position = 0; position < bwt_size; ++position) {
    if (position != eob_position) {
      uint32 ch = static_cast<byte>(bwt[position]);
      uint32 rank = count[ch + 1];
      uint32 rank_low24 = rank & 0x00FFFFFF;
      //uint32 rank_hi8 = rank >> 24;
      bwt_rank_low24[position] = (ch << 24) + rank_low24;
      if (0 == rank_low24) {
        rank_milestone_buffer[ch].push_back(position);
      }
      ++count[ch + 1];
    }
  }
  std::partial_sum(count.begin(), count.end(), count.begin());
  assert(count[256] == bwt_size);

  for (int ch = 0; ch < 256; ++ch) {
    rank_milestone_buffer[ch].push_back(bwt_size);
  }

  std::vector<byte>* result = AllocateMemory(bwt_size - 1);
  uint32 index = 0;

  uint32 position = 0;
  while (position != eob_position) {
    uint32 bwt_and_rank = bwt_rank_low24[position];
    uint32 ch = bwt_and_rank >> 24;
    (*result)[index++] = ch;
    uint32 rank_low24 = bwt_and_rank & 0x00FFFFFF;
    uint32 rank_hi8 = 0;
    while (rank_milestone_buffer[ch][rank_hi8 + 1] <= position) ++rank_hi8;
    position = count[ch] + (rank_hi8 << 24) + rank_low24;
  }
  // If the BWT or the EOB position contain errors (or are garbage),
  // it is likely that the cycle ends prematurely.
  // In this case, we return false.
  //assert(index == bwt_size);
  return result;
}

} //namespace bwtc
