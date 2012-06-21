/**
 * @file GrammarTest.cpp
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
 * Tests for grammar.
 */


#include "../preprocessors/Grammar.hpp"
#include "TestStreams.hpp"

#define BOOST_TEST_MODULE 
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <string>
#include <vector>

namespace bwtc {
int verbosity = 0;

namespace tests {

BOOST_AUTO_TEST_SUITE(GrammarSerialization)

BOOST_AUTO_TEST_CASE(EmptyGrammar) {
  TestStream stream;
  Grammar original;
  original.writeGrammar(&stream);
  stream.reset();
  Grammar copy;
  copy.readGrammar(&stream);
  BOOST_CHECK_EQUAL(copy.numberOfRules(), original.numberOfRules());
  BOOST_CHECK_EQUAL(copy.numberOfSpecialSymbols(),
                    original.numberOfSpecialSymbols());
  BOOST_CHECK_EQUAL(copy.numberOfRules(), 0);
  BOOST_CHECK_EQUAL(copy.numberOfSpecialSymbols(), 0);
}

BOOST_AUTO_TEST_CASE(SpecialSymbols) {
  TestStream stream;
  Grammar original;
  original.beginUpdatingRules();

  std::vector<uint16> specialPairs;
  std::vector<byte> freed, spec;
  freed.push_back('a'); freed.push_back('b');
  spec.push_back('X'); spec.push_back('Y');

  original.addRule('a', '1','2');
  original.expandAlphabet(freed, spec, specialPairs);
  original.endUpdatingRules();

  BOOST_CHECK_EQUAL((specialPairs[0] >> 8) & 0xff, 'X');
  BOOST_CHECK_EQUAL(specialPairs[0] & 0xff, 'Y');
  BOOST_CHECK_EQUAL((specialPairs[1] >> 8) & 0xff, 'Y');
  BOOST_CHECK_EQUAL(specialPairs[1] & 0xff, 'X');
  
  BOOST_CHECK_EQUAL(original.numberOfRules(), 1);
  BOOST_CHECK_EQUAL(original.numberOfSpecialSymbols(), 2);

  original.writeGrammar(&stream);
  stream.reset();
  Grammar copy;
  copy.readGrammar(&stream);

  BOOST_CHECK_EQUAL(copy.numberOfRules(), original.numberOfRules());
  BOOST_CHECK_EQUAL(copy.numberOfSpecialSymbols(),
                    original.numberOfSpecialSymbols());
}

BOOST_AUTO_TEST_SUITE_END()

} //namespace tests
} //namespace bwtc

