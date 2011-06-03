/**
 * @file long_sequences.cc
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
 * Header for detecting and replacing long repetitive sequences.
 */

#include "../globaldefs.hpp"
#include "longsequences.h"
#include "sequence_detector.h"
#include <vector>

namespace bwtc {
namespace long_sequences {

typedef hash_functions::MaskHasher Hasher;
const uint32 LONG_SEQUENCE = 1000;
const uint32 MIN_PERIOD = 10;
const uint32 error_val = Max<uint32>::max;

typedef chunk bucket_struct;

void FormBucketsAndNamePositions(std::vector<uint32> *names,
                                 std::vector<chunk> *chunks,
                                 std::vector<bucket_struct> *buckets);

uint32 CompressPositionOrdered(const std::vector<uint32>& bucket_starts,
                               std::vector<chunk> *pos_ordered);

void SortIntoBuckets(std::vector<uint32> *hash_values,
                     std::vector<chunk> *chunks,
                     std::vector<bucket_struct> *buckets);

uint32 GetPosition(const bucket_struct *b);

void SetBucketBeginFlag(bucket_struct *b);

bool StartOfBucket(const bucket_struct *b);

void StringSort(bucket_struct *begin, bucket_struct *end,
                const std::vector<chunk>& pos_ordered,
                uint32 pre_len, byte *from, uint32 str_len);

uint32 NextValidPos(const std::vector<uint32>& bucket_starts, uint32 pos);

void SortBuckets(byte *from, uint32 win_length,
                 const std::vector<uint32>& bucket_starts,
                 std::vector<bucket_struct> *buckets,
                 const std::vector<chunk>& pos_ordered);

void HandleBuckets(byte *from, uint32 win_size,
                   std::vector<bucket_struct> *buckets,
                   std::vector<chunk> *pos_ordered);


uint64 CompressSequences(byte *from, uint32 length, uint32 window_size);


} //namespace long_sequences
} //namespace bwtc
