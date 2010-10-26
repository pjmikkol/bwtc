/**
 * @file expander_test.cc
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
 * Unit tests for bwtc::long_sequences::Expander-class. Sorting and TODO WHAT?
 * 
 */

#define BOOST_TEST_MODULE 
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include <cassert>
#include <iostream>
#include <vector>
#include "../globaldefs.h"
#include "../preprocessors/longsequences.h"

#define BWTC_SEQ_DET_SCAN_CONSTS_
namespace bwtc {
namespace long_sequences {
const uint32 kMinPeriod = 4;
const uint32 kWindowSize = 8;
}
}

#include "../preprocessors/expander.h"

namespace bwtc {
int verbosity = 0;

namespace tests {

using namespace long_sequences;

/**
 * Vectors needed for bucket sort
 */
struct Vectors {
 public:
  Vectors() { }
  std::vector<uint32> uints;
  std::vector<chunk> chunks;
  std::vector<chunk> buckets;
};

BOOST_AUTO_TEST_SUITE(BucketSort)

BOOST_AUTO_TEST_CASE(BucketSort) {
  Vectors v;
  v.uints.reserve(6);
  v.uints.push_back(2); v.uints.push_back(0); v.uints.push_back(2);
  v.uints.push_back(1); v.uints.push_back(4); v.uints.push_back(1);
  v.chunks.reserve(10);
  v.chunks.push_back(chunk(1,0));
  v.chunks.push_back(chunk(14,0));
  v.chunks.push_back(chunk(20,2));
  v.chunks.push_back(chunk(25,5));
  v.chunks.push_back(chunk(30,4));
  v.chunks.push_back(chunk(35,4));
  v.chunks.push_back(chunk(40,4));
  v.chunks.push_back(chunk(45,3));
  v.chunks.push_back(chunk(50,4));
  v.chunks.push_back(chunk(58,2));
  std::vector<uint32> *b = SortIntoBuckets(&v.uints, &v.chunks, &v.buckets);
  BOOST_CHECK_EQUAL(b->size(), 3); BOOST_CHECK_EQUAL((*b)[0], 2);
  BOOST_CHECK_EQUAL((*b)[1], 4); BOOST_CHECK_EQUAL((*b)[2], 8); 
  BOOST_CHECK_EQUAL(v.buckets[0].position, 14);
  BOOST_CHECK_EQUAL(v.buckets[1].position, 1);
  BOOST_CHECK_EQUAL(v.buckets[2].position, 58);
  BOOST_CHECK_EQUAL(v.buckets[3].position, 20);
  BOOST_CHECK_EQUAL(v.buckets[4].position, 50);
  BOOST_CHECK_EQUAL(v.buckets[5].position, 40);
  BOOST_CHECK_EQUAL(v.buckets[6].position, 35);
  BOOST_CHECK_EQUAL(v.buckets[7].position, 30);
  delete b;
}

BOOST_AUTO_TEST_SUITE_END()

} //namespace tests
} //namespace long_sequences
