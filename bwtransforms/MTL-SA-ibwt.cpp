/**
 * @file MTL-SA-ibwt.hpp
 * @author Juha Karkkainen <juha.karkkainen@cs.helsinki.fi>
 * @author Dominik Kempa <dominik.kempa@cs.helsinki.fi>
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
 * Implementation of the MTL-SA algorithm for inverting BWT.
 */

#include <cassert>
#include <numeric>  // For partial_sum.

#include "../globaldefs.hpp"
#include "MTL-SA-ibwt.hpp"

namespace bwtc {

uint64 MTL_SA_InverseBWTransform::maxBlockSize(uint64 memory_budget) const {
  return memory_budget / sizeof(uint32);
}

std::vector<byte>* MTL_SA_InverseBWTransform::doTransform(
    const byte* bwt, uint64 bwt_size, uint64 eob_position) {
  // Original string cannot be empty.
  assert(bwt_size >= 2);

  // We use 'data' to store:
  //   a) LF^2[i] (4 bytes),
  //   b) a pair (bwt[LF[i]], bwt[i]) (2 bytes)
  // for each position i = 0, .., bwt_size - 1. In total the array takes roughly
  // 6 * bwt_size bytes. Having this array enables restoring the original string
  // two characters at a time.
  // Information a) and b) for position i is always accessed together. To reduce
  // the number of cache misses the following layout is used:
  //
  //   layout  | - - - - | - - | - - | - - - - | - - - - | - - | - - | - - - -|
  //   bytes        4       2     2       4         4       2     2       4
  //   meaning   LF^2[0]   P[0]  P[1]  LF^2[1]   LF^2[2]   P[2]  P[3]  LF^2[3]
  //
  // Where P[i] is a pair (bwt[LF[i]], bwt[i]).
  std::vector<uint32> data(3 * ((bwt_size + 1) / 2));

  // pairs_count[code] counts the number of occurrences of pair AB in the
  // original string, where code is the byte sequence: A, B interpreted as 16
  // bit integer.
  std::vector<uint32> pairs_count(256 * 256 + 1, 0);

  // EOB_pairs stores the numbers of rows from BWT matrix that ends with a pair
  // containing the EOB symbol.
  std::pair<uint32, uint32> EOB_pairs = 
    computePairs(bwt, bwt_size, eob_position, data, pairs_count);

  // Fill the remaining information in the 'data' array.
  computeLF2(bwt_size, EOB_pairs, pairs_count, data);

  std::vector<byte>* result = allocateMemory(bwt_size - 1);
  uint16 *result_ptr = (uint16 *)(&(*result)[0]);
  uint32 filled = 0;
  uint32 *data_ptr = &data[0];
  uint32 position = 0;

  // Number of full pairs to restore.
  uint32 aligned_end = (bwt_size - 1) / 2;

  // Restore the original string, two characters at a time.
  while (filled < aligned_end) {
    // base and offset are separated to save time during address computation.
    uint32 base = 3 * (position / 2);
    uint32 offset = position & 1;
    uint32 next_position = data[base + (offset << 1)];
    uint16 char_pair = ((uint16 *)(data_ptr + base + 1))[offset];
    result_ptr[filled++] = char_pair;
    position = next_position;
  }

  // Single character left.
  if (!(bwt_size & 1)) {
    uint32 base = 3 * (position / 2);
    uint32 offset = position & 1;
    uint16 char_pair = ((uint16 *)(data_ptr + base + 1))[offset];
    (*result)[bwt_size - 2] = ((byte *)&char_pair)[0];
  }

  return result;
}

std::pair<uint32, uint32> MTL_SA_InverseBWTransform::computePairs(
    const byte* bwt, uint64 bwt_size, uint64 eob_position,
    std::vector<uint32> &data,  std::vector<uint32> &pairs_count) {
  // During the algorithm we create an array of bytes which is later
  // interpreted as an array of 16 bit integers, hence we need to know
  // the machine endianness to correctly encode the byte pairs.
  uint16 value = 1;
  bool little_endian = ((unsigned char *)&value)[0];

  // Used to keep track of characters in the first column on BWT matrix.
  std::vector<uint32> count(256 + 1, 0);

  // Count EOB.
  count[0] = 1;

  // Count other characters.
  for (uint32 position = 0; position < bwt_size; ++position) {
    if (position != eob_position) {
      uint32 ch = bwt[position];
      ++count[ch + 1];
    }
  }
  std::partial_sum(count.begin(), count.end(), count.begin());
  assert(count[256] == bwt_size);

  // Used to compute LF mapping during sequential scan of BWT.
  std::vector<uint32> rank(256 + 1);
  std::copy(count.begin(), count.end(), rank.begin());

  // During the scan of BWT stores the corresponding character from the
  // first column of BWT matrix.
  int32 first_column_ch = 0;
 
  uint16 *pairs_ptr = (uint16 *)(&data[0] + 1);
  uint32 first_EOB_pair = eob_position;
  uint32 second_EOB_pair = 0;

  // Handle position 0.
  uint32 bwt_current = bwt[0];
  uint32 LF_current = rank[bwt_current];
  uint32 bwt_previous = bwt[LF_current];
  ++rank[bwt_current];
  if (little_endian) {
    *pairs_ptr = bwt_current + (bwt_previous << 8);
  } else {
    *pairs_ptr = (bwt_current << 8) + bwt_previous;
  }
  ++pairs_ptr;

  // Handle positions [1, bwt_size).
  if (little_endian) {
    for (uint32 position = 1; position < bwt_size; ++position) {
      if (position != eob_position) {
        bwt_current = bwt[position];
        LF_current = rank[bwt_current];
        if (LF_current == eob_position) {
          second_EOB_pair = position;
        }
        bwt_previous = bwt[LF_current];
        ++rank[bwt_current];
        *pairs_ptr = bwt_current + (bwt_previous << 8);
        while (count[first_column_ch + 1] <= position) {
          ++first_column_ch;
        }
        ++pairs_count[(first_column_ch << 8) + bwt_current + 1];
      }
      pairs_ptr += 1 + ((position & 1) << 2);
    }
    // Transpose pairs_count.
    for (uint32 low = 0; low < 256; ++low) {
      for (uint32 high = low + 1; high < 256; ++high) {
        std::swap(pairs_count[(low << 8) + high + 1],
                  pairs_count[(high << 8) + low + 1]);
      }
    }
  } else {
    for (uint32 position = 1; position < bwt_size; ++position) {
      if (position != eob_position) {
        bwt_current = bwt[position];
        LF_current = rank[bwt_current];
        if (LF_current == eob_position) {
          second_EOB_pair = position;
        }
        bwt_previous = bwt[LF_current];
        ++rank[bwt_current];      
        *pairs_ptr = (bwt_current << 8) + bwt_previous;
        while (count[first_column_ch + 1] <= position) {
          ++first_column_ch;
        }
        ++pairs_count[(first_column_ch << 8) + bwt_current + 1];
      }
      pairs_ptr += 1 + ((position & 1) << 2);
    }
  }

  // Count the only pair starting with an EOB.
  ++pairs_count[0];

  // Count the only pair ending with an EOB.
  bwt_current = bwt[0];
  if (little_endian) {
    ++pairs_count[bwt_current << 8];
  } else {
    ++pairs_count[(255 << 8) + bwt_current];
  }

  return std::make_pair(first_EOB_pair, second_EOB_pair);
}

void MTL_SA_InverseBWTransform::computeLF2(uint64 bwt_size,
    std::pair<uint32, uint32> &EOB_pairs, std::vector<uint32> &pairs_count,
    std::vector<uint32> &data) {
  // Accumulate the counts, so we can use it to compute LF^2. This is analogous
  // to MTL algorithm, except now we use the counts for pairs rather than single
  // characters.
  std::partial_sum(pairs_count.begin(), pairs_count.end(), pairs_count.begin());

  uint32 *data_ptr = &data[0];
  uint32 bwt_index = 0;
  uint32 *data_end = data_ptr + 3 * (bwt_size / 2);
  uint16 *pairs_ptr = (uint16 *)(&data[0] + 1);

  // Two entries per one loop execution are computed. This is due to nonlinear
  // layout in which powers of LF mapping are stored.
  while (data_ptr < data_end) {
    *data_ptr = pairs_count[*pairs_ptr];
    if (bwt_index != EOB_pairs.first &&
        bwt_index != EOB_pairs.second) {
      ++pairs_count[*pairs_ptr];
    }
    data_ptr += 2;
    ++pairs_ptr;
    ++bwt_index;
    *data_ptr = pairs_count[*pairs_ptr];
    if (bwt_index != EOB_pairs.first &&
        bwt_index != EOB_pairs.second) {
      ++pairs_count[*pairs_ptr];
    }
    ++data_ptr;
    pairs_ptr += 5;
    ++bwt_index;
  }
  if (bwt_size & 1) {
    *data_ptr = pairs_count[*pairs_ptr]++;
  }
}

} //namespace bwtc

