/**
 * @file UtilsTest.cpp
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
 * Unit tests for various functionality located at preprocessor.cc.
 */

#define BOOST_TEST_MODULE 
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include "../globaldefs.hpp"

#include "../Utils.hpp"

namespace bwtc {
int verbosity = 0;

namespace tests {

//TODO: test also other functions
using namespace utils;

BOOST_AUTO_TEST_SUITE(GlobalFunctions)

BOOST_AUTO_TEST_CASE(MostSigBit) {
  BOOST_CHECK_EQUAL(0x80, mostSignificantBit(0xff));
  BOOST_CHECK_EQUAL(0x00, mostSignificantBit(0x00));
  BOOST_CHECK_EQUAL(0x400000, mostSignificantBit(0x500032));
  BOOST_CHECK_EQUAL(0x1000000, mostSignificantBit(0x14f3da2));
  BOOST_CHECK_EQUAL(0x200, mostSignificantBit(0x3f3));
}

BOOST_AUTO_TEST_CASE(MostSigBit16) {
  BOOST_CHECK_EQUAL(0x8000, mostSignificantBit16(0xffff));
  BOOST_CHECK_EQUAL(0x00, mostSignificantBit16(0x00));
  BOOST_CHECK_EQUAL(0x4000, mostSignificantBit16(0x5032));
  BOOST_CHECK_EQUAL(0x1000, mostSignificantBit16(0x14f3));
  BOOST_CHECK_EQUAL(0x200, mostSignificantBit16(0x3f3));
}

BOOST_AUTO_TEST_CASE(LogFlr) {
  BOOST_CHECK_EQUAL(7, logFloor(0xff));
  BOOST_CHECK_EQUAL(26, logFloor(0x61adf3f));
  BOOST_CHECK_EQUAL(30, logFloor(0x561adf3f));
  BOOST_CHECK_EQUAL(31, logFloor(0xf61adf3f));
}

BOOST_AUTO_TEST_CASE(RunFrequencyCounts1) {
  const char * str = "aabacca";
  uint64 freqs[256] = {0};
  calculateRunFrequencies(freqs, (const byte*)str, strlen(str));
  BOOST_CHECK_EQUAL(freqs['a'], 3);
  BOOST_CHECK_EQUAL(freqs['c'], 1);
  BOOST_CHECK_EQUAL(freqs['b'], 1);
}

BOOST_AUTO_TEST_CASE(RunFrequencyCounts2) {
  const char * str = "wwwaaaawawawabbbb";
  uint64 freqs[256] = {0};
  calculateRunFrequencies(freqs, (const byte*)str, strlen(str));
  BOOST_CHECK_EQUAL(freqs['w'], 4);
  BOOST_CHECK_EQUAL(freqs['a'], 4);
  BOOST_CHECK_EQUAL(freqs['b'], 1);
}

BOOST_AUTO_TEST_SUITE_END()

} //namespace tests
} //namespace long_sequences

