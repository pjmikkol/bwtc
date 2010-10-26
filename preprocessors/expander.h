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

std::vector<uint32> *SortIntoBuckets(std::vector<uint32> *hash_values,
                                    std::vector<chunk> *chunks,
                                    std::vector<chunk> *buckets);
void FormBucketsAndNamePositions(std::vector<uint32> *names,
                                 std::vector<chunk> *chunks,
                                 std::vector<chunk> *buckets);

} //namespace long_sequences
} //namespace bwtc


#endif
