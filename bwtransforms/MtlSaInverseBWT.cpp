/**
 * @file MtlSaInverseBWT.hpp
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
#include "MtlSaInverseBWT.hpp"

namespace bwtc {

uint64 MtlSaInverseBWTransform::maxBlockSize(uint64 memory_budget) const {
  return memory_budget / sizeof(uint32);
}

void computeData(const byte *bwt, uint64 bwt_size, uint32 *data,
    uint32 eob_position) {
  // The function assumes that bwt contains the EOB symbol, which conceptually
  // is the smallest character in the alphabet but never takes any action that
  // depends on its value, hence actual value of EOB is irrelevant, in
  // particular it can be > 0.

  // During the algorithm we create an array of bytes which is later
  // interpreted as an array of 16 bit integers, hence we need to know
  // the machine endianness to correctly encode the byte pairs.
  uint16 value = 1;
  bool little_endian = ((unsigned char *)&value)[0];

  // Used to keep track of characters in the first column on BWT matrix.
  std::vector<uint32> count(256 + 1, 0);
  uint32 *count_ptr = &count[0];
  
  // Used to compute LF mapping during sequential scan of BWT.
  std::vector<uint32> rank(256 + 1);
  uint32 *rank_ptr = &rank[0];
  
  // Used to compute the values of LF^2.
  std::vector<uint32> pairs_count(256 * 256 + 1, 0);
  uint32 *pairs_count_ptr = &pairs_count[0];

  // Count EOB.
  count_ptr[0] = 1;

  // Count other characters.
  for (uint32 position = 0; position < bwt_size; ++position) {
    if (position != eob_position) {
      uint32 ch = bwt[position];
      ++count_ptr[ch + 1];
    }
  }
  std::partial_sum(count.begin(), count.end(), count.begin());
  std::copy(count.begin(), count.end(), rank.begin());
  assert(count[256] == bwt_size);

  // During the scan of BWT stores the corresponding character from the
  // first column of BWT matrix.
  int32 first_column_ch = 0;
  uint16 *pairs_ptr = (uint16 *)(&data[0] + 1);

  // Store the indices of BWT matrix rows that contain the EOB symbol in the
  // last 2 columns. These positions usually must be handled separately.
  uint32 first_eob_pair = eob_position;
  uint32 second_eob_pair = 0;

  // Handle position 0.
  //   a) compute previous character
  uint32 bwt_current = bwt[0];
  uint32 LF_0 = rank_ptr[bwt_current];
  uint32 bwt_previous = bwt[LF_0];
  //   b) store pair
  if (little_endian) {
    *pairs_ptr = bwt_current + (bwt_previous << 8);
  } else {
    *pairs_ptr = (bwt_current << 8) + bwt_previous;
  }
  ++pairs_ptr;
  //   c) update rank
  ++rank_ptr[bwt_current];

  // Handle positions [1, bwt_size).
  if (little_endian) {
    for (uint32 position = 1; position < eob_position; ++position) {
      //   a) compute previous character
      bwt_current = bwt[position];
      uint32 LF_current = rank_ptr[bwt_current];
      bwt_previous = bwt[LF_current];
      //   b) store pair
      *pairs_ptr = bwt_current + (bwt_previous << 8);
      pairs_ptr += 1 + ((position & 1) << 2);
      //   c) update rank
      ++rank_ptr[bwt_current];
      if (LF_current == eob_position) {
        second_eob_pair = position;
      }
      while (count_ptr[first_column_ch + 1] <= position) {
        ++first_column_ch;
      }
      ++pairs_count_ptr[(first_column_ch << 8) + bwt_current + 1];
    }
    *pairs_ptr = (bwt[0] << 8);
    pairs_ptr += 1 + ((eob_position & 1) << 2);
    for (uint32 position = eob_position + 1; position < bwt_size; ++position) {
      bwt_current = bwt[position];
      uint32 LF_current = rank_ptr[bwt_current];
      bwt_previous = bwt[LF_current];
      *pairs_ptr = bwt_current + (bwt_previous << 8);
      pairs_ptr += 1 + ((position & 1) << 2);
      ++rank_ptr[bwt_current];
      if (LF_current == eob_position) {
        second_eob_pair = position;
      }
      while (count_ptr[first_column_ch + 1] <= position) {
        ++first_column_ch;
      }
      ++pairs_count_ptr[(first_column_ch << 8) + bwt_current + 1];
    }
    // Transpose pairs_count.
    for (uint32 low = 0; low < 256; ++low) {
      for (uint32 high = low + 1; high < 256; ++high) {
        std::swap(pairs_count[(low << 8) + high + 1],
                  pairs_count[(high << 8) + low + 1]);
      }
    }
  } else {
    for (uint32 position = 1; position < eob_position; ++position) {
      bwt_current = bwt[position];
      uint32 LF_current = rank_ptr[bwt_current];
      bwt_previous = bwt[LF_current];
      *pairs_ptr = (bwt_current << 8) + bwt_previous;
      pairs_ptr += 1 + ((position & 1) << 2);
      ++rank_ptr[bwt_current];
      if (LF_current == eob_position) {
        second_eob_pair = position;
      }
      while (count_ptr[first_column_ch + 1] <= position) {
        ++first_column_ch;
      }
      ++pairs_count_ptr[(first_column_ch << 8) + bwt_current + 1];
    }
    *pairs_ptr = bwt[0];
    pairs_ptr += 1 + ((eob_position & 1) << 2);
    for (uint32 position = eob_position + 1; position < bwt_size; ++position) {
      bwt_current = bwt[position];
      uint32 LF_current = rank_ptr[bwt_current];
      bwt_previous = bwt[LF_current];
      *pairs_ptr = (bwt_current << 8) + bwt_previous;
      pairs_ptr += 1 + ((position & 1) << 2);
      ++rank_ptr[bwt_current];      
      if (LF_current == eob_position) {
        second_eob_pair = position;
      }
      while (count_ptr[first_column_ch + 1] <= position) {
        ++first_column_ch;
      }
      ++pairs_count_ptr[(first_column_ch << 8) + bwt_current + 1];
    }
  }

  // Count the only pair starting with an EOB.
  ++pairs_count_ptr[0];

  // Count the only pair ending with an EOB.
  bwt_current = bwt[0];
  if (little_endian) {
    ++pairs_count_ptr[bwt_current << 8];
  } else {
    ++pairs_count_ptr[(255 << 8) + bwt_current];
  }

  // Accumulate the counts, so we can use it to compute LF^2.
  std::partial_sum(pairs_count.begin(), pairs_count.end(), pairs_count.begin());

  uint32 *data_ptr = &data[0];
  uint32 bwt_index = 0;
  uint32 *data_end = data_ptr + 3 * (bwt_size / 2);
  pairs_ptr = (uint16 *)(&data[0] + 1);

  // Two entries per one loop execution are computed. This is due to nonlinear
  // layout in which powers of LF are stored.
  while (data_ptr < data_end) {
    // 1
    if (bwt_index == first_eob_pair) {
      *data_ptr = LF_0;
    } else if (bwt_index == second_eob_pair) {
      *data_ptr = 0;
    } else {
      *data_ptr = pairs_count_ptr[*pairs_ptr];
      ++pairs_count_ptr[*pairs_ptr];
    }
    data_ptr += 2;
    ++pairs_ptr;
    ++bwt_index;
    // 2
    if (bwt_index == first_eob_pair) {
      *data_ptr = LF_0;
    } else if (bwt_index == second_eob_pair) {
      *data_ptr = 0;
    } else {
      *data_ptr = pairs_count_ptr[*pairs_ptr];
      ++pairs_count_ptr[*pairs_ptr];
    }
    ++data_ptr;
    pairs_ptr += 5;
    ++bwt_index;
  }
  
  // Handle the remaining single position.
  if (bwt_size & 1) {
    if (bwt_index == first_eob_pair) {
      *data_ptr = LF_0;
    } else if (bwt_index == second_eob_pair) {
      *data_ptr = 0;
    } else {
      *data_ptr = pairs_count_ptr[*pairs_ptr];
      ++pairs_count_ptr[*pairs_ptr];
    }
  }
}

std::vector<byte>* MtlSaInverseBWTransform::doTransform(const byte* bwt,
    uint64 bwt_size, const std::vector<uint32> &LFpowers) {
  assert(bwt_size >= 2);
  assert(LFpowers.size() > 0);
  uint32 eob_position = LFpowers[0];

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
  uint32 *data = new uint32[3 * ((bwt_size + 1) / 2)];
  computeData(bwt, bwt_size, data, eob_position);

  std::vector<byte> *result = allocateMemory(bwt_size - 1);
  uint32 *data_ptr = &data[0];
  byte *result_ptr = &(*result)[0];
  uint32 starting_positions = LFpowers.size();
  uint32 block_size = bwt_size / starting_positions;
  uint32 to_restore = bwt_size - 1;

  // Stores the set of current LF powers (one per block).
  uint32 positions[starting_positions];
  std::copy(LFpowers.begin(), LFpowers.end(), positions);

  // If the block size is small, don't deploy parallel inversion.
  if (block_size <= 1) {
    starting_positions = 1;
    block_size = bwt_size;
  }

  // Stores pointers to positions in text that are about to be restored.
  uint16 *dest_ptr[starting_positions];
  for (uint32 i = 0; i < starting_positions; ++i) {
    dest_ptr[i] = (uint16 *)(result_ptr + i * block_size - 1);
  }

  // Restore the first pair from each block.
  for (uint32 block_id = 0; block_id < starting_positions; ++block_id) {
    uint32 position = positions[block_id];
    uint32 base = 3 * (position / 2);
    uint32 offset = position & 1;
    uint32 next_position = data_ptr[base + (offset << 1)];
    uint16 char_pair = ((uint16 *)(data_ptr + base + 1))[offset];
    positions[block_id] = next_position;
    if (block_id == 0) {
      // Skip the EOB symbol.
      ++dest_ptr[block_id];
      result_ptr[0] = ((unsigned char *)&char_pair)[1];
    } else { 
      // Not the first pair, decode as normal.
      *dest_ptr[block_id]++ = char_pair;
    }
  }

  // Restore the main part of each block, two characters at a time,
  // simultaneously from multiple starting positions.
  for (uint32 filled = 1; filled < block_size / 2; ++filled) {
    for (uint32 block_id = 0; block_id < starting_positions; ++block_id) {
      uint32 position = positions[block_id];
      uint32 base = 3 * (position / 2);
      uint32 offset = position & 1;
      uint32 next_position = data_ptr[base + (offset << 1)];
      uint16 char_pair = ((uint16 *)(data_ptr + base + 1))[offset];
      positions[block_id] = next_position;
      *dest_ptr[block_id]++ = char_pair;
    }
  }

  // If the block size is odd then there is one remaining character in each
  // block. This loop takes case of that. The last block is handled separately
  // because it might contain more than one remaining character.
  if (block_size & 1) {
    for (uint32 block_id = 0; block_id + 1 < starting_positions; ++block_id) {
      uint32 position = positions[block_id];
      uint32 base = 3 * (position / 2);
      uint32 offset = position & 1;
      uint16 char_pair = ((uint16 *)(data_ptr + base + 1))[offset];
      result_ptr[(block_id + 1) * block_size - 2] =
        ((unsigned char *)&char_pair)[0];
    }
  }

  // Restore the remaining characters from the last (longest) block,
  // possibly except the last character, if the last block has odd length.
  uint32 index = (starting_positions - 1) * block_size - 1 + 2 * (block_size / 2);
  while (index + 1 < to_restore) {
    uint32 position = positions[starting_positions - 1];
    uint32 base = 3 * (position / 2);
    uint32 offset = position & 1;
    uint32 next_position = data_ptr[base + (offset << 1)];
    uint16 char_pair = ((uint16 *)(data_ptr + base + 1))[offset];
    positions[starting_positions - 1] = next_position;
    *dest_ptr[starting_positions - 1]++ = char_pair;
    index += 2;
  }

  // Single character left in the last block.
  if (index != to_restore) {
    uint32 position = positions[starting_positions - 1];
    uint32 base = 3 * (position / 2);
    uint32 offset = position & 1;
    uint16 char_pair = ((uint16 *)(data_ptr + base + 1))[offset];
    result_ptr[to_restore - 1] = ((unsigned char *)&char_pair)[0];
  }

  delete[] data;
  return result;
}

} //namespace bwtc

