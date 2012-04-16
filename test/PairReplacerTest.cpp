/**
 * @file PairReplacerTest.cpp
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
 * Tests for PairReplacer.
 */


#include "../preprocessors/PairReplacer.hpp"

#define BOOST_TEST_MODULE 
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <string>
#include <vector>

namespace bwtc {
int verbosity = 0;

namespace tests {

BOOST_AUTO_TEST_SUITE(analyseTests)

BOOST_AUTO_TEST_CASE(analyseAtOnce) {
  std::string abab = "ab";
  std::vector<byte> data;
  for(size_t j = 0; j < 10000; ++j) {
    for(size_t i = 0; i < abab.size(); ++i) {
      data.push_back(abab[i]);
    }
  }
  for(size_t i = 0; i < 256; ++i) {
    if (i == 'a' || i == 'b') continue;
    data.push_back((byte) i);
  }
  PairReplacer pr(true);
  pr.analyseData(&data[0], data.size());
  size_t rep = pr.decideReplacements();
  BOOST_CHECK_EQUAL(rep, 1);
  // headerSize == 5, totalSize == 10000+254+5
  std::vector<byte> result;
  result.resize(10259);
  size_t hSize = pr.writeHeader(&result[0]);
  BOOST_CHECK_EQUAL(hSize,5);
  size_t cSize = pr.writeReplacedVersion(&data[0], data.size(), &result[5]);
  BOOST_CHECK_EQUAL(cSize,10254);
}

BOOST_AUTO_TEST_SUITE_END()


} //namespace tests
} //namespace bwtc

