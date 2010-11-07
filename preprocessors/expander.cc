/**
 * @file expander.cc
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

#include "../globaldefs.h"
#include "longsequences.h"
#include "expander.h"
#include <vector>
#include <utility>

namespace bwtc {
namespace long_sequences {


std::vector<uint32> *SortIntoBuckets(std::vector<uint32> *hash_values,
                                     std::vector<chunk> *chunks,
                                     std::vector<chunk> *buckets)
{
  uint32 total = 0;
  std::vector<uint32> *bucket_limits = new std::vector<uint32>();
  for(std::vector<uint32>::iterator it = hash_values->begin();
      it != hash_values->end(); ++it)
  {
    if(*it > 1) {
      total += *it; *it = total;
      bucket_limits->push_back(total);
    }
    else *it = Max<uint32>::max;
  }
  buckets->resize(total);
  FormBucketsAndNamePositions(hash_values, chunks, buckets);
  return bucket_limits;
}

void FormBucketsAndNamePositions(std::vector<uint32> *names,
                                 std::vector<chunk> *chunks,
                                 std::vector<chunk> *buckets) 
{
  uint32 position = 0;
  std::vector<uint32>& n_vec = *names;
  std::vector<chunk>& buckets_ = *buckets;
  for(std::vector<chunk>::iterator it=chunks->begin(); it != chunks->end();
      ++it)
  {
    if(n_vec[it->hash_value] != Max<uint32>::max) 
      buckets_[--n_vec[it->hash_value]] = chunk(it->position, position++);
  }
}

} //namespace bwtc
} //namespace long_sequences
