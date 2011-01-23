/**
 * @file preprocessor_tests.cc
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

#include "../globaldefs.h"

#include "../utils.h"

namespace bwtc {
int verbosity = 0;

namespace tests {

//TODO: test also other functions
using namespace utils;

BOOST_AUTO_TEST_SUITE(GlobalFunctions)

BOOST_AUTO_TEST_CASE(MostSigBit) {
  BOOST_CHECK_EQUAL(0x80, MostSignificantBit(0xff));
  BOOST_CHECK_EQUAL(0x00, MostSignificantBit(0x00));
  BOOST_CHECK_EQUAL(0x400000, MostSignificantBit(0x500032));
  BOOST_CHECK_EQUAL(0x1000000, MostSignificantBit(0x14f3da2));
  BOOST_CHECK_EQUAL(0x200, MostSignificantBit(0x3f3));
}

BOOST_AUTO_TEST_CASE(MostSigBit16) {
  BOOST_CHECK_EQUAL(0x8000, MostSignificantBit16(0xffff));
  BOOST_CHECK_EQUAL(0x00, MostSignificantBit16(0x00));
  BOOST_CHECK_EQUAL(0x4000, MostSignificantBit16(0x5032));
  BOOST_CHECK_EQUAL(0x1000, MostSignificantBit16(0x14f3));
  BOOST_CHECK_EQUAL(0x200, MostSignificantBit16(0x3f3));
}

BOOST_AUTO_TEST_CASE(LogFlr) {
  BOOST_CHECK_EQUAL(7, LogFloor(0xff));
  BOOST_CHECK_EQUAL(26, LogFloor(0x61adf3f));
  BOOST_CHECK_EQUAL(30, LogFloor(0x561adf3f));
  BOOST_CHECK_EQUAL(31, LogFloor(0xf61adf3f));
}

BOOST_AUTO_TEST_SUITE_END()

} //namespace tests
} //namespace long_sequences

