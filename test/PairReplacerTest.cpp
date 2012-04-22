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
  result.resize(10254+9);
  size_t cSize = pr.writeReplacedVersion(&data[0], data.size(), &result[0]);
  BOOST_CHECK_EQUAL(cSize,10254);
  uint32 gSize = grammar.writeGrammar(&result[10254]);
  BOOST_CHECK_EQUAL(grammar.numberOfSpecialSymbols(),0);
  BOOST_CHECK_EQUAL(grammar.numberOfRules(),1);
  BOOST_CHECK_EQUAL(gSize,9);

  PostProcessor post(true);
  post.postProcess(&result);
  BOOST_CHECK_EQUAL(result.size(), data.size());
  for(size_t i = 0; i < result.size(); ++i)
    BOOST_CHECK_EQUAL(result[i], data[i]);
}

BOOST_AUTO_TEST_CASE(pairOfSame) {
  std::cout << std::endl << std::endl << std::endl;
  //Should replace aa with a without freeing symbols
  Grammar grammar;
  std::string aa = "aa";
  std::vector<byte> data;
  data.resize(20000+255);
  uint32 k = 0;
  for(size_t j = 0; j < 10000; ++j) {
    for(size_t i = 0; i < aa.size(); ++i) {
      data[k++] = (byte)aa[i];
    }
  }
  for(size_t i = 0; i < 256; ++i) {
    if (i == 'a') continue;
    data[k++] = ((byte) i);
  }
  PairReplacer pr(grammar, true);
  pr.analyseData(&data[0], data.size());
  size_t rep = pr.decideReplacements();
  BOOST_CHECK_EQUAL(rep, 1);
  // size of replaced string == 20000 + 255
  std::vector<byte> result;
  result.resize(10000 + 255 + 9);
  size_t cSize = pr.writeReplacedVersion(&data[0], data.size(), &result[0]);
  BOOST_CHECK_EQUAL(cSize,10255);
  uint32 gSize = grammar.writeGrammar(&result[10255]);
  BOOST_CHECK_EQUAL(grammar.numberOfSpecialSymbols(),0);
  BOOST_CHECK_EQUAL(grammar.numberOfRules(),1);
  BOOST_CHECK_EQUAL(gSize,9);

  grammar.printRules();
  
  PostProcessor post(true);
  post.postProcess(&result);
  BOOST_CHECK_EQUAL(result.size(), data.size());
  for(size_t i = 0; i < result.size(); ++i)
    BOOST_CHECK_EQUAL(result[i], data[i]);
}

BOOST_AUTO_TEST_CASE(twoForFree) {
  std::cout << "twoForFree" << std::endl;
  // c and d are used as variables after freeing pair cd
  Grammar grammar;
  std::string aa = "acdf";
  std::vector<byte> data;
  data.resize(40000+256);
  uint32 k = 0;
  for(size_t j = 0; j < 10000; ++j) {
    for(size_t i = 0; i < aa.size(); ++i) {
      data[k++] = (byte)aa[i];
    }
  }
  for(size_t i = 0; i < 256; ++i) {
    data[k++] = ((byte) i);
  }
  PairReplacer pr(grammar, true);
  pr.analyseData(&data[0], data.size());
  size_t rep = pr.decideReplacements();
  BOOST_CHECK_EQUAL(rep, 2);

  std::vector<byte> result;
  result.resize(20001 + 255 + 13);
  size_t cSize = pr.writeReplacedVersion(&data[0], data.size(), &result[0]);
  BOOST_CHECK_EQUAL(cSize,20256);
  uint32 gSize = grammar.writeGrammar(&result[20256]);
  BOOST_CHECK_EQUAL(grammar.numberOfSpecialSymbols(),0);
  BOOST_CHECK_EQUAL(grammar.numberOfRules(),2);
  BOOST_CHECK_EQUAL(gSize,13);

  grammar.printRules();
  
  PostProcessor post(true);
  post.postProcess(&result);
  BOOST_CHECK_EQUAL(result.size(), data.size());
  for(size_t i = 0; i < result.size(); ++i)
    BOOST_CHECK_EQUAL(result[i], data[i]);
}


BOOST_AUTO_TEST_CASE(simpleSpecialScenario) {
  std::cout << "\n\nSimple special scenario" << std::endl;
  Grammar grammar;
  std::string aa = "aceh";
  std::vector<byte> data;
  data.resize(40000+256);
  uint32 k = 0;
  for(size_t j = 0; j < 10000; ++j) {
    for(size_t i = 0; i < aa.size(); ++i) {
      data[k++] = (byte)aa[i];
    }
  }
  for(size_t i = 0; i < 256; ++i) {
    data[k++] = ((byte) i);
  }
  PairReplacer pr(grammar, true);
  pr.analyseData(&data[0], data.size());
  size_t rep = pr.decideReplacements();
  BOOST_CHECK_EQUAL(rep, 2);
  // size of replaced string == 20000 + 256 + 4 (two specials, two freed) 
  std::vector<byte> result;
  result.resize(20000 + 256 + 4 + 17);
  size_t cSize = pr.writeReplacedVersion(&data[0], data.size(), &result[0]);
  BOOST_CHECK_EQUAL(cSize,20260);
  uint32 gSize = grammar.writeGrammar(&result[20260]);
  BOOST_CHECK_EQUAL(grammar.numberOfSpecialSymbols(),2);
  BOOST_CHECK_EQUAL(grammar.numberOfRules(),2);
  BOOST_CHECK_EQUAL(gSize,17);

  grammar.printRules();
  
  PostProcessor post(true);
  post.postProcess(&result);
  BOOST_CHECK_EQUAL(result.size(), data.size());
  for(size_t i = 0; i < result.size(); ++i)
    BOOST_CHECK_EQUAL(result[i], data[i]);
}


BOOST_AUTO_TEST_CASE(freeVariableFromReplacements) {
  std::cout << std::endl << std::endl << std::endl;
  //Note that pair replacer doesn't that a could be used for free after two
  //replacements of pairs
  Grammar grammar;
  std::string ab = "acad";
  std::vector<byte> data;
  data.resize(40000+255);
  uint32 k = 0;
  for(size_t j = 0; j < 10000; ++j) {
    for(size_t i = 0; i < ab.size(); ++i) {
      data[k++] = (byte)ab[i];
    }
  }
  for(size_t i = 0; i < 256; ++i) {
    if (i == 'a') continue;
    data[k++] = ((byte) i);
  }
  PairReplacer pr(grammar, true);
  pr.analyseData(&data[0], data.size());
  size_t rep = pr.decideReplacements();
  BOOST_CHECK_EQUAL(rep, 2);
  // size of replaced string == 20000 + 254
  std::vector<byte> result;
  result.resize(20000 + 258 + 16);
  size_t cSize = pr.writeReplacedVersion(&data[0], data.size(), &result[0]);
  BOOST_CHECK_EQUAL(cSize,20258);
  uint32 gSize = grammar.writeGrammar(&result[20258]);
  BOOST_CHECK_EQUAL(grammar.numberOfSpecialSymbols(),2);
  BOOST_CHECK_EQUAL(grammar.numberOfRules(),2);
  BOOST_CHECK_EQUAL(gSize,16);

  grammar.printRules();
  
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
  result.resize(10259+ 12);
  size_t cSize = pr.writeReplacedVersion(&data[0], data.size(), &result[0]);
  BOOST_CHECK_EQUAL(cSize,10259);
  uint32 gSize = grammar.writeGrammar(&result[10259]);
  BOOST_CHECK_EQUAL(grammar.numberOfSpecialSymbols(),2);
  BOOST_CHECK_EQUAL(grammar.numberOfRules(),1);
  BOOST_CHECK_EQUAL(gSize,12);

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
  result.resize(10512+ 13);
  size_t cSize = pr.writeReplacedVersion(&data[0], data.size(), &result[0]);
  BOOST_CHECK_EQUAL(cSize,10512);
  uint32 gSize = grammar.writeGrammar(&result[10512]);
  BOOST_CHECK_EQUAL(grammar.numberOfSpecialSymbols(),2);
  BOOST_CHECK_EQUAL(grammar.numberOfRules(),1);
  BOOST_CHECK_EQUAL(gSize,13);
  BOOST_CHECK_EQUAL(result.back(),1);

  PostProcessor post(true);
  post.postProcess(&result);
  BOOST_CHECK_EQUAL(result.size(), data.size());
  for(size_t i = 0; i < result.size(); ++i) 
    BOOST_CHECK_EQUAL(result[i], data[i]);

}


BOOST_AUTO_TEST_CASE(MoreReplacements) {
  std::cout << std::endl << std::endl << std::endl;

  Grammar grammar;
  std::string ab = "aereaeab";
  std::vector<byte> data;
  data.resize(80000+510);
  uint32 k = 0;
  for(size_t j = 0; j < 10000; ++j) {
    for(size_t i = 0; i < ab.size(); ++i) {
      data[k++] = (byte)ab[i];
    }
  }
  for(size_t j = 0; j < 2; ++j) {
    for(size_t i = 0; i < 256; ++i) {
      if(i == 'a') continue;
      data[k++] = ((byte) i);
    }
  }
  BOOST_CHECK_EQUAL(k, data.size());
  PairReplacer pr(grammar, true);
  pr.analyseData(&data[0], data.size());
  size_t rep = pr.decideReplacements();
  BOOST_CHECK_EQUAL(rep, 3);
  std::vector<byte> result;
  result.resize(80000);
  size_t cSize = pr.writeReplacedVersion(&data[0], data.size(), &result[0]);
  uint32 gSize = grammar.writeGrammar(&result[cSize]);
  BOOST_CHECK_EQUAL(grammar.numberOfRules(),3);
  BOOST_CHECK_EQUAL(grammar.numberOfSpecialSymbols(),2);
  result.resize(cSize+gSize);

  grammar.printRules();
  
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

BOOST_AUTO_TEST_CASE(NewSpecialOnRightSide) {
  std::cout << std::endl << std::endl << std::endl;

  Grammar grammar;
  std::string ab = "af";
  std::vector<byte> data;
  data.resize(20000+512);
  uint32 k = 0;
  for(size_t j = 0; j < 10000; ++j) {
    for(size_t i = 0; i < ab.size(); ++i) {
      data[k++] = (byte)ab[i];
    }
  }
  for(size_t j = 0; j < 2; ++j) {
    for(size_t i = 0; i < 256; ++i) {
      data[k++] = ((byte) i);
    }
  }
  BOOST_CHECK_EQUAL(k, data.size());
  PairReplacer pr(grammar, true);
  pr.analyseData(&data[0], data.size());
  size_t rep = pr.decideReplacements();
  BOOST_CHECK_EQUAL(rep, 1);
  std::vector<byte> result;
  result.resize(10000 + 512 + 6 + 12);
  size_t cSize = pr.writeReplacedVersion(&data[0], data.size(), &result[0]);
  uint32 gSize = grammar.writeGrammar(&result[cSize]);
  BOOST_CHECK_EQUAL(grammar.numberOfRules(),1);
  BOOST_CHECK_EQUAL(grammar.numberOfSpecialSymbols(),2);
  BOOST_CHECK_EQUAL(cSize, 10518);
  BOOST_CHECK_EQUAL(gSize, 12);
  result.resize(10530);

  grammar.printRules();
  
  PostProcessor post(true);
  post.postProcess(&result);
  BOOST_CHECK_EQUAL(result.size(), data.size());
  for(size_t i = 0; i < result.size(); ++i) 
    BOOST_CHECK_EQUAL(result[i], data[i]);

}

BOOST_AUTO_TEST_CASE(twoRounds) {
  std::cout << std::endl << std::endl << std::endl;

  Grammar grammar;
  std::string ab = "afad";
  std::vector<byte> data;
  data.resize(40000+511);
  uint32 k = 0;
  for(size_t j = 0; j < 10000; ++j) {
    for(size_t i = 0; i < ab.size(); ++i) {
      data[k++] = (byte)ab[i];
    }
  }
  for(size_t j = 0; j < 2; ++j) {
    for(size_t i = 0; i < 256; ++i) {
      if(i == 'a') continue;
      data[k++] = ((byte) i);
    }
  }
  data[k++] = 'a';
  BOOST_CHECK_EQUAL(k, data.size());

  std::vector<byte> result[2];
  size_t cSize;
  for(size_t i = 0; i < 2; ++i) {
    PairReplacer pr(grammar, true);
    std::vector<byte> *v;
    if(i == 0) v = &data;
    else v = &result[i-1];
    
    pr.analyseData(&(*v)[0], v->size());

    size_t rep = pr.decideReplacements();

    result[i].resize(40000);
    cSize = pr.writeReplacedVersion(&(*v)[0], v->size(), &result[i][0]);
    result[i].resize(cSize);
  }
  result[1].resize(result[1].size() + 1000);
  uint32 gSize = grammar.writeGrammar(&result[1][cSize]);
  result[1].resize(cSize+gSize);
  grammar.printRules();
  
  PostProcessor post(true);
  post.postProcess(&result[1]);
  BOOST_CHECK_EQUAL(result[1].size(), data.size());
  for(size_t i = 0; i < result[1].size(); ++i) 
    BOOST_CHECK_EQUAL(result[1][i], data[i]);


}

BOOST_AUTO_TEST_SUITE_END()


} //namespace tests
} //namespace bwtc

