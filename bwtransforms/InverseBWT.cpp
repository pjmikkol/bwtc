/**
 * @file InverseBWT.cpp
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
 * Implementation for the inverse of BWT. This particular implementation
 * is released originally here:
 * @see http://code.google.com/p/dcs-bwt-compressor/
 */

#include <cassert>

#include <numeric>  // for partial_sum
#include <vector>

#include "../globaldefs.hpp"
#include "InverseBWT.hpp"
#include "MtlSaInverseBWT.hpp"
#include "../BWTBlock.hpp"
#include "../Profiling.hpp"

namespace bwtc {

InverseBWTransform* giveInverseTransformer() {
  //return new FastInverseBWTransform();
  return new MtlSaInverseBWTransform();
}

void InverseBWTransform::doTransform(BWTBlock& block) {
  byte *data = block.begin();
  *block.end() = data[block.LFpowers()[0]];
  doTransform(block.begin(), block.size()+1, block.LFpowers());
}

uint64 FastInverseBWTransform::maxBlockSize(uint64 memory_budget) const {
  memory_budget -= kMemoryOverhead;
  return memory_budget / sizeof(uint32);
}

void FastInverseBWTransform::doTransform(
    byte* bwt, uint32 bwt_size, const std::vector<uint32>& LFpowers)
{
  PROFILE("FastInverseBWTransform::doTransform");
  uint32 eob_position = LFpowers[0];
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

  uint32 index = 0;

  uint32 position = 0;
  while (position != eob_position) {
    uint32 bwt_and_rank = bwt_rank_low24[position];
    uint32 ch = bwt_and_rank >> 24;
    bwt[index++] = ch;
    uint32 rank_low24 = bwt_and_rank & 0x00FFFFFF;
    uint32 rank_hi8 = 0;
    while (rank_milestone_buffer[ch][rank_hi8 + 1] <= position) ++rank_hi8;
    position = count[ch] + (rank_hi8 << 24) + rank_low24;
  }
  // If the BWT or the EOB position contain errors (or are garbage),
  // it is likely that the cycle ends prematurely.
  // In this case, we return false.
  //assert(index == bwt_size);
}

} //namespace bwtc
