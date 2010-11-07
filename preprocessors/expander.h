/**
 * @file expander.h
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
 */

#ifndef BWTC_EXPANDER_H_
#define BWTC_EXPANDER_H_

#include "../globaldefs.h"
#include "longsequences.h"
#include <vector>

namespace bwtc {
namespace long_sequences {


/** Bucket-sorts hash values and gives new names to positions.
 *
 * @param hash_values An array with counts of occurences of each hash value.
 *                    The counts are indexed by their hash value.
 * @param chunks Vector which consists of chunks found at the scanning phase.
 * @param buckets The result is storded to here. This is assumed to be empty
 *                when called this function.
 * @return An array which holds the information of the buckets. i'th
 *         hash values is stored to range [a[i-1], a[i]) where a is the
 *         vector returned.
 */
std::vector<uint32> *SortIntoBuckets(std::vector<uint32> *hash_values,
                                    std::vector<chunk> *chunks,
                                    std::vector<chunk> *buckets);

void FormBucketsAndNamePositions(std::vector<uint32> *names,
                                 std::vector<chunk> *chunks,
                                 std::vector<chunk> *buckets);

/** Sorts bucket of chunks at the range [begin, end).
 *
 * @param begin Start of the range to be sorted.
 * @param end End of the range to be sorted.
 * @param sub_buckets Stores the found sub buckets in this vector. Expected to
 *                    be empty when called.
 * @param source Since chunks are always associated to some source we need
 *               this for conducting the comparisons. 
 */
void SortBucket(chunk *begin, chunk *end, std::vector<uint32> *sub_buckets,
                byte *source);

void HandleSubBucket(chunk *begin, chunk *end, byte *source);


} //namespace long_sequences
} //namespace bwtc


#endif
