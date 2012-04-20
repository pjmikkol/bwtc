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
#include "../preprocessors/Grammar.hpp"
#include "../preprocessors/Postprocessor.hpp"

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
  Grammar grammar;
  std::string ab = "ab";
  std::vector<byte> data;
  data.resize(20000+254);
  uint32 k = 0;
  for(size_t j = 0; j < 10000; ++j) {
    for(size_t i = 0; i < ab.size(); ++i) {
      data[k++] = (byte)ab[i];
    }
  }
  for(size_t i = 0; i < 256; ++i) {
    if (i == 'a' || i == 'b') continue;
    data[k++] = ((byte) i);
  }
  PairReplacer pr(grammar, true);
  pr.analyseData(&data[0], data.size());
  size_t rep = pr.decideReplacements();
  BOOST_CHECK_EQUAL(rep, 1);
  // size of replaced string == 10000+254
  std::vector<byte> result;
  result.resize(10254+8);
  size_t cSize = pr.writeReplacedVersion(&data[0], data.size(), &result[0]);
  BOOST_CHECK_EQUAL(cSize,10254);
  uint32 gSize = grammar.writeGrammar(&result[10254]);
  BOOST_CHECK_EQUAL(grammar.numberOfSpecialSymbols(),0);
  BOOST_CHECK_EQUAL(grammar.numberOfRules(),1);
  BOOST_CHECK_EQUAL(gSize,8);


  PostProcessor post(true);
  post.postProcess(&result);
  BOOST_CHECK_EQUAL(result.size(), data.size());
  for(size_t i = 0; i < result.size(); ++i)
    BOOST_CHECK_EQUAL(result[i], data[i]);
}

BOOST_AUTO_TEST_CASE(MakeSymbolsFree) {
  std::cout << std::endl << std::endl << std::endl;
  Grammar grammar;
  std::string ab = "ac";
  std::vector<byte> data;
  data.resize(20000+256);
  uint32 k = 0;
  for(size_t j = 0; j < 10000; ++j) {
    for(size_t i = 0; i < ab.size(); ++i) {
      data[k++] = (byte)ab[i];
    }
  }
  for(size_t i = 0; i < 256; ++i) {
    data[k++] = ((byte) i);
  }
  PairReplacer pr(grammar, true);
  pr.analyseData(&data[0], data.size());
  size_t rep = pr.decideReplacements();
  BOOST_CHECK_EQUAL(rep, 1);
  // size of replaced string == 10000+253 + 2 + 2 + 2
  std::vector<byte> result;
  //grammar header = 1 + 3 + 1 + 1 + 2 + 1+2
  result.resize(10259+ 11);
  size_t cSize = pr.writeReplacedVersion(&data[0], data.size(), &result[0]);
  BOOST_CHECK_EQUAL(cSize,10259);
  uint32 gSize = grammar.writeGrammar(&result[10259]);
  BOOST_CHECK_EQUAL(grammar.numberOfSpecialSymbols(),2);
  BOOST_CHECK_EQUAL(grammar.numberOfRules(),1);
  BOOST_CHECK_EQUAL(gSize,11);

  PostProcessor post(true);
  post.postProcess(&result);
  BOOST_CHECK_EQUAL(result.size(), data.size());
  for(size_t i = 0; i < result.size(); ++i) 
    BOOST_CHECK_EQUAL(result[i], data[i]);
}

BOOST_AUTO_TEST_CASE(MakeVariableFree) {
  std::cout << std::endl << std::endl << std::endl;

  Grammar grammar;
  std::string ab = "ac";
  std::vector<byte> data;
  data.resize(20000+509);
  uint32 k = 0;
  for(size_t j = 0; j < 10000; ++j) {
    for(size_t i = 0; i < ab.size(); ++i) {
      data[k++] = (byte)ab[i];
    }
  }
  for(size_t j = 0; j < 2; ++j) {
    for(size_t i = 0; i < 256; ++i) {
      if(i == 'X' || i == 'Y' || i == 'a') continue;
      data[k++] = ((byte) i);
    }
  }
  data[k++] = 'a';
  data[k++] = 'X';
  data[k++] = 'Y';
  BOOST_CHECK_EQUAL(k, data.size());
  PairReplacer pr(grammar, true);
  pr.analyseData(&data[0], data.size());
  size_t rep = pr.decideReplacements();
  BOOST_CHECK_EQUAL(rep, 1);
  // size of replaced string == 10000+506 + 3*2
  std::vector<byte> result;
  //grammar header = 1 + 1 + 2 + 2 + 1+1+2
  result.resize(10512+ 12);
  size_t cSize = pr.writeReplacedVersion(&data[0], data.size(), &result[0]);
  BOOST_CHECK_EQUAL(cSize,10512);
  uint32 gSize = grammar.writeGrammar(&result[10512]);
  BOOST_CHECK_EQUAL(grammar.numberOfSpecialSymbols(),2);
  BOOST_CHECK_EQUAL(grammar.numberOfRules(),1);
  BOOST_CHECK_EQUAL(gSize,12);
  BOOST_CHECK_EQUAL(result.back(),1);

  PostProcessor post(true);
  post.postProcess(&result);
  BOOST_CHECK_EQUAL(result.size(), data.size());
  for(size_t i = 0; i < result.size(); ++i) 
    BOOST_CHECK_EQUAL(result[i], data[i]);

}

BOOST_AUTO_TEST_CASE(MultipleRounds) {
  std::cout << std::endl << std::endl << std::endl;

  Grammar grammar;
  std::string ab = "ac";
  std::vector<byte> data;
  data.resize(20000+509);
  uint32 k = 0;
  for(size_t j = 0; j < 10000; ++j) {
    for(size_t i = 0; i < ab.size(); ++i) {
      data[k++] = (byte)ab[i];
    }
  }
  for(size_t j = 0; j < 2; ++j) {
    for(size_t i = 0; i < 256; ++i) {
      if(i == 'X' || i == 'Y' || i == 'a') continue;
      data[k++] = ((byte) i);
    }
  }
  data[k++] = 'a';
  data[k++] = 'X';
  data[k++] = 'Y';
  PairReplacer pr(grammar, true);
  pr.analyseData(&data[0], data.size());
  size_t rep = pr.decideReplacements();
  BOOST_CHECK_EQUAL(rep, 1);
  // size of replaced string == 10000+506 + 3*2
  std::vector<byte> result;
  //grammar header = 1 + 1 + 2 + 2 + 1+1+2
  result.resize(10512+ 11);
  size_t cSize = pr.writeReplacedVersion(&data[0], data.size(), &result[0]);
  BOOST_CHECK_EQUAL(cSize,10512);

  PairReplacer pr1(grammar, true);
  pr1.analyseData(&result[0], 10512);
  rep = pr1.decideReplacements();
  BOOST_CHECK_EQUAL(rep, 1);
  // size of replaced string == 5000+ 505 + 4*2
  //grammar header = 1 + 1 + 2 + 3 + 1+2+4
  BOOST_CHECK_EQUAL(grammar.isSpecial(result[0]), false);
  std::vector<byte> res2;
  res2.resize(5512 + 40);
  cSize = pr1.writeReplacedVersion(&result[0], 10512, &res2[0]);
  BOOST_CHECK_EQUAL(cSize,5512);
  

  uint32 gSize = grammar.writeGrammar(&res2[5512]);
  BOOST_CHECK_EQUAL(grammar.numberOfSpecialSymbols(),2);
  BOOST_CHECK_EQUAL(grammar.numberOfRules(),2);
  BOOST_CHECK_EQUAL(res2[cSize+gSize-1],2);

  res2.resize(gSize + cSize);

  grammar.printRules();

  PostProcessor post(true);
  post.postProcess(&res2);
  BOOST_CHECK_EQUAL(res2.size(), data.size());
  for(size_t i = 0; i < res2.size(); ++i) 
    BOOST_CHECK_EQUAL(res2[i], data[i]);

}

BOOST_AUTO_TEST_SUITE_END()


} //namespace tests
} //namespace bwtc

