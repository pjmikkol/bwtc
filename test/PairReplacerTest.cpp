/**
 * @file PairReplacerTest.cpp
 * @author Pekka Mikkola <pmikkol@gmail.com>
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

BOOST_AUTO_TEST_SUITE(analyseTests)

BOOST_AUTO_TEST_CASE(analyseAtOnce) {
  std::cout << "analyseAtOnce" << std::endl;
  std::vector<byte> data;
  data.resize(20255);
  std::string ab = "ac";
  uint32 k = 0;
  for(size_t j = 0; j < 10000; ++j) {
    for(size_t i = 0; i < ab.size(); ++i) {
      data[k++] = (byte)ab[i];
    }
  }
  for(size_t i = 0; i < 256; ++i) {
    if (i == 'g') continue;
    data[k++] = (byte) i;
  }

  Grammar grammar;
  PairReplacer pr(grammar, true);
  pr.analyseData(&data[0], data.size());
  size_t rep = pr.decideReplacements();
  BOOST_CHECK_EQUAL(rep, 1);

  std::vector<byte> result(data);

  std::vector<byte> vec;
  TestStream serializedGrammar(vec);
  grammar.writeGrammar(&serializedGrammar);
  serializedGrammar.reset();
  size_t cSize = pr.writeReplacedVersion(&result[0], result.size());

  // size of replaced string == 10000+255
  BOOST_CHECK_EQUAL(cSize,10255);
  BOOST_CHECK_EQUAL(grammar.numberOfSpecialSymbols(),0);
  BOOST_CHECK_EQUAL(grammar.numberOfRules(),1);

  
  Grammar postGrammar;
  postGrammar.readGrammar(&serializedGrammar);
  for(size_t i = 0; i < 256; ++i)
    BOOST_CHECK_EQUAL(grammar.frequencies()[i],postGrammar.frequencies()[i]);

  Postprocessor post(true, postGrammar);

  std::vector<byte> original;
  post.uncompress(&result[0], cSize, original);

  BOOST_CHECK_EQUAL(original.size(), data.size());
  for(size_t i = 0; i < original.size(); ++i) {
    BOOST_CHECK_EQUAL(original[i], data[i]);

  }
}

BOOST_AUTO_TEST_CASE(pairOfSame) {
  std::cout << "pairOfSame" << std::endl;
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
    data[k++] = (byte) i;
  }

  Grammar grammar;
  PairReplacer pr(grammar, true);
  pr.analyseData(&data[0], data.size());
  size_t rep = pr.decideReplacements();
  BOOST_CHECK_EQUAL(rep, 1);

  std::vector<byte> result(data);
  std::vector<byte> vec;
  TestStream serializedGrammar(vec);
  grammar.writeGrammar(&serializedGrammar);
  serializedGrammar.reset();
  size_t cSize = pr.writeReplacedVersion(&result[0], result.size());

  // size of replaced string == 10000 + 255 + 3
  BOOST_CHECK_EQUAL(cSize,10258);
  BOOST_CHECK_EQUAL(grammar.numberOfSpecialSymbols(),2);
  BOOST_CHECK_EQUAL(grammar.numberOfRules(),1);

  Grammar postGrammar;
  postGrammar.readGrammar(&serializedGrammar);
  for(size_t i = 0; i < 256; ++i)
    BOOST_CHECK_EQUAL(grammar.frequencies()[i],postGrammar.frequencies()[i]);

  Postprocessor post(true, postGrammar);

  std::vector<byte> original;
  post.uncompress(&result[0], cSize, original);


  //BOOST_CHECK_EQUAL(original.size(), data.size());
  for(size_t i = 0; i < original.size(); ++i) {
    BOOST_CHECK_EQUAL(original[i], data[i]);
  }
}

BOOST_AUTO_TEST_CASE(twoForFree) {
  std::cout << "twoForFree" << std::endl;
  Grammar grammar;
  std::string aa = "acdf";
  std::vector<byte> data;
  data.resize(40256);
  uint32 k = 0;
  for(size_t j = 0; j < 10000; ++j) {
    for(size_t i = 0; i < aa.size(); ++i) {
      data[k++] = (byte)aa[i];
    }
  }
  for(size_t i = 0; i < 256; ++i) {
    data[k++] = (byte) i;
  }
  PairReplacer pr(grammar, true);
  pr.analyseData(&data[0], data.size());
  size_t rep = pr.decideReplacements();
  BOOST_CHECK_EQUAL(rep, 2);

  std::vector<byte> result(data);

  std::vector<byte> vec;
  TestStream serializedGrammar(vec);
  grammar.writeGrammar(&serializedGrammar);
  serializedGrammar.reset();
  size_t cSize = pr.writeReplacedVersion(&result[0], result.size());

  // Size of result == 20256 + 2*2
  BOOST_CHECK_EQUAL(cSize,20260);
  BOOST_CHECK_EQUAL(grammar.numberOfSpecialSymbols(),2);
  BOOST_CHECK_EQUAL(grammar.numberOfRules(),2);

  Grammar postGrammar;
  postGrammar.readGrammar(&serializedGrammar);
  Postprocessor post(true, postGrammar);

  std::vector<byte> original;
  post.uncompress(&result[0], cSize, original);

  BOOST_CHECK_EQUAL(original.size(), data.size());
  for(size_t i = 0; i < original.size(); ++i)
    BOOST_CHECK_EQUAL(original[i], data[i]);
}


BOOST_AUTO_TEST_CASE(simpleSpecialScenario) {
  std::cout << "simpleSpecialScenario" << std::endl;
  std::string aa = "aceh";
  std::vector<byte> data;
  data.resize(40256);
  uint32 k = 0;
  for(size_t j = 0; j < 10000; ++j) {
    for(size_t i = 0; i < aa.size(); ++i) {
      data[k++] = (byte)aa[i];
    }
  }
  for(size_t i = 0; i < 256; ++i) {
    data[k++] = ((byte) i);
  }
  
  Grammar grammar;
  PairReplacer pr(grammar, true);
  pr.analyseData(&data[0], data.size());
  size_t rep = pr.decideReplacements();
  BOOST_CHECK_EQUAL(rep, 2);

  // size of replaced string == 20000 + 256 + 4 (two specials, two freed) 
  std::vector<byte> result(data);

  std::vector<byte> vec;
  TestStream serializedGrammar(vec);
  grammar.writeGrammar(&serializedGrammar);
  serializedGrammar.reset();
  size_t cSize = pr.writeReplacedVersion(&result[0], result.size());

  BOOST_CHECK_EQUAL(cSize,20260);
  BOOST_CHECK_EQUAL(grammar.numberOfSpecialSymbols(),2);
  BOOST_CHECK_EQUAL(grammar.numberOfRules(),2);

  Grammar postGrammar;
  postGrammar.readGrammar(&serializedGrammar);
  Postprocessor post(true, postGrammar);

  std::vector<byte> original;
  post.uncompress(&result[0], cSize, original);

  BOOST_CHECK_EQUAL(original.size(), data.size());
  for(size_t i = 0; i < original.size(); ++i)
    BOOST_CHECK_EQUAL(original[i], data[i]);
}


BOOST_AUTO_TEST_CASE(freeCharacterFromReplacements) {
  std::cout << "freeCharacterFromReplacements" << std::endl;
  std::string ab = "acad";
  std::vector<byte> data;
  data.resize(40255);
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

  Grammar grammar;
  PairReplacer pr(grammar, true);
  pr.analyseData(&data[0], data.size());
  size_t rep = pr.decideReplacements();
  BOOST_CHECK_EQUAL(rep, 2);

  // size of replaced string == 20000 + 255 + 4
  std::vector<byte> result(data);

  std::vector<byte> vec;
  TestStream serializedGrammar(vec);
  grammar.writeGrammar(&serializedGrammar);
  serializedGrammar.reset();
  size_t cSize = pr.writeReplacedVersion(&result[0], result.size());

  BOOST_CHECK_EQUAL(cSize,20259);
  BOOST_CHECK_EQUAL(grammar.numberOfSpecialSymbols(),2);
  BOOST_CHECK_EQUAL(grammar.numberOfRules(),2);

  Grammar postGrammar;
  postGrammar.readGrammar(&serializedGrammar);
  Postprocessor post(true, postGrammar);

  std::vector<byte> original;
  post.uncompress(&result[0], cSize, original);

  BOOST_CHECK_EQUAL(original.size(), data.size());
  for(size_t i = 0; i < original.size(); ++i)
    BOOST_CHECK_EQUAL(original[i], data[i]);
}

BOOST_AUTO_TEST_CASE(MakeSymbolsFree) {
  std::cout << "MakeSymbolsFree" << std::endl;
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
    data[k++] = (byte) i;
  }

  Grammar grammar;
  PairReplacer pr(grammar, true);
  pr.analyseData(&data[0], data.size());
  size_t rep = pr.decideReplacements();
  BOOST_CHECK_EQUAL(rep, 1);

  // size of replaced string == 10000 + 253 + 3x2
  std::vector<byte> result(data);

  std::vector<byte> vec;
  TestStream serializedGrammar(vec);
  grammar.writeGrammar(&serializedGrammar);
  serializedGrammar.reset();
  size_t cSize = pr.writeReplacedVersion(&result[0], result.size());

  BOOST_CHECK_EQUAL(cSize,10259);
  BOOST_CHECK_EQUAL(grammar.numberOfSpecialSymbols(),2);
  BOOST_CHECK_EQUAL(grammar.numberOfRules(),1);

  Grammar postGrammar;
  postGrammar.readGrammar(&serializedGrammar);
  Postprocessor post(true, postGrammar);

  std::vector<byte> original;
  post.uncompress(&result[0], cSize, original);

  BOOST_CHECK_EQUAL(original.size(), data.size());
  for(size_t i = 0; i < original.size(); ++i)
    BOOST_CHECK_EQUAL(original[i], data[i]);
}

BOOST_AUTO_TEST_CASE(MakeVariableFree) {
  std::cout << "MakeVariableFree" << std::endl;
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

  Grammar grammar;
  PairReplacer pr(grammar, true);
  pr.analyseData(&data[0], data.size());
  size_t rep = pr.decideReplacements();
  BOOST_CHECK_EQUAL(rep, 1);

  // size of replaced string == 10000+506 + 2 + 1 + 2*2
  std::vector<byte> result(data);

  std::vector<byte> vec;
  TestStream serializedGrammar(vec);
  grammar.writeGrammar(&serializedGrammar);
  serializedGrammar.reset();
  size_t cSize = pr.writeReplacedVersion(&result[0], result.size());

  BOOST_CHECK_EQUAL(cSize,10513);
  BOOST_CHECK_EQUAL(grammar.numberOfSpecialSymbols(),2);
  BOOST_CHECK_EQUAL(grammar.numberOfRules(),1);

  Grammar postGrammar;
  postGrammar.readGrammar(&serializedGrammar);
  Postprocessor post(true, postGrammar);

  std::vector<byte> original;
  post.uncompress(&result[0], cSize, original);

  BOOST_CHECK_EQUAL(original.size(), data.size());
  for(size_t i = 0; i < original.size(); ++i)
    BOOST_CHECK_EQUAL(original[i], data[i]);
}


BOOST_AUTO_TEST_CASE(MoreReplacements) {
  std::cout << "MoreReplacements" << std::endl;
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

  Grammar grammar;
  PairReplacer pr(grammar, true);
  pr.analyseData(&data[0], data.size());
  size_t rep = pr.decideReplacements();
  BOOST_CHECK_EQUAL(rep, 3);

  std::vector<byte> result(data);

  std::vector<byte> vec;
  TestStream serializedGrammar(vec);
  grammar.writeGrammar(&serializedGrammar);
  serializedGrammar.reset();

  size_t cSize = pr.writeReplacedVersion(&result[0], result.size());

  BOOST_CHECK_EQUAL(grammar.numberOfRules(),3);
  BOOST_CHECK_EQUAL(grammar.numberOfSpecialSymbols(),3);

  Grammar postGrammar;
  postGrammar.readGrammar(&serializedGrammar);
  Postprocessor post(true, postGrammar);

  std::vector<byte> original;
  post.uncompress(&result[0], cSize, original);

  BOOST_CHECK_EQUAL(original.size(), data.size());
  for(size_t i = 0; i < original.size(); ++i)
    BOOST_CHECK_EQUAL(original[i], data[i]);
}

BOOST_AUTO_TEST_CASE(MultipleRounds) {
  std::cout << "MultipleRounds" << std::endl;
  std::string ab = "ac";
  std::vector<byte> data;
  data.resize(20000+508);
  uint32 k = 0;
  for(size_t j = 0; j < 10000; ++j) {
    for(size_t i = 0; i < ab.size(); ++i) {
      data[k++] = (byte)ab[i];
    }
  }
  for(size_t j = 0; j < 2; ++j) {
    for(size_t i = 0; i < 256; ++i) {
      if(i == 'X' || i == 'Y' || i == 'a') continue;
      data[k++] = (byte) i;
    }
  }
  data[k++] = 'X';
  data[k++] = 'Y';

  Grammar grammar;
  PairReplacer pr(grammar, true);
  pr.analyseData(&data[0], data.size());
  size_t rep = pr.decideReplacements();
  BOOST_CHECK_EQUAL(rep, 1);

  // size of replaced string == 10000+504 + 2x2 + 2x2
  std::vector<byte> result(data);

  size_t cSize = pr.writeReplacedVersion(&result[0], result.size());

  BOOST_CHECK_EQUAL(cSize,10512);
  PairReplacer pr1(grammar, true);
  pr1.analyseData(&result[0], 10512);
  rep = pr1.decideReplacements();
  BOOST_CHECK_EQUAL(rep, 1);

  // Should use 'a' as a new variable
  // size of replaced string == 5000 + 504 + 2x2 + 2x2
  BOOST_CHECK_EQUAL(grammar.isSpecial(result[0]), false);
  cSize = pr1.writeReplacedVersion(&result[0], cSize);

  std::vector<byte> vec;
  TestStream serializedGrammar(vec);
  grammar.writeGrammar(&serializedGrammar);
  serializedGrammar.reset();

  
  BOOST_CHECK(cSize == 5512);
  
  
  BOOST_CHECK_EQUAL(grammar.numberOfSpecialSymbols(),2);
  BOOST_CHECK_EQUAL(grammar.numberOfRules(),2);

  Grammar postGrammar;
  postGrammar.readGrammar(&serializedGrammar);
  Postprocessor post(true, postGrammar);

  std::vector<byte> original;
  post.uncompress(&result[0], cSize, original);

  BOOST_CHECK_EQUAL(original.size(), data.size());
  for(size_t i = 0; i < original.size(); ++i)
    BOOST_CHECK_EQUAL(original[i], data[i]);
}

BOOST_AUTO_TEST_CASE(NewSpecialOnRightSide) {
  std::cout << "NewSpecialOnRightSide" << std::endl;
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
      data[k++] = (byte) i;
    }
  }
  BOOST_CHECK_EQUAL(k, data.size());

  PairReplacer pr(grammar, true);
  pr.analyseData(&data[0], data.size());
  size_t rep = pr.decideReplacements();
  BOOST_CHECK_EQUAL(rep, 1);

  std::vector<byte> result(data);
  std::vector<byte> vec;
  TestStream serializedGrammar(vec);
  grammar.writeGrammar(&serializedGrammar);
  serializedGrammar.reset();
  size_t cSize = pr.writeReplacedVersion(&result[0], result.size());

  BOOST_CHECK_EQUAL(grammar.numberOfRules(),1);
  BOOST_CHECK_EQUAL(grammar.numberOfSpecialSymbols(),2);
  BOOST_CHECK_EQUAL(cSize, 10000 + 512 + 2*3);

  Grammar postGrammar;
  postGrammar.readGrammar(&serializedGrammar);
  Postprocessor post(true, postGrammar);

  std::vector<byte> original;
  post.uncompress(&result[0], cSize, original);

  BOOST_CHECK_EQUAL(original.size(), data.size());
  for(size_t i = 0; i < original.size(); ++i)
    BOOST_CHECK_EQUAL(original[i], data[i]);
}

BOOST_AUTO_TEST_CASE(twoRounds) {
  std::cout << "twoRounds" << std::endl;
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

  std::vector<byte> orig(data);
  size_t cSize = data.size();
  for(size_t i = 0; i < 2; ++i) {
    PairReplacer pr(grammar, true);
    pr.analyseData(&data[0], cSize);

    pr.decideReplacements();

    cSize = pr.writeReplacedVersion(&data[0], cSize);
  }
  std::vector<byte> vec;
  TestStream serializedGrammar(vec);
  grammar.writeGrammar(&serializedGrammar);
  serializedGrammar.reset();

  Grammar postGrammar;
  postGrammar.readGrammar(&serializedGrammar);
  Postprocessor post(true, postGrammar);

  std::vector<byte> original;
  post.uncompress(&data[0], cSize, original);

  BOOST_CHECK_EQUAL(original.size(), orig.size());
  for(size_t i = 0; i < original.size(); ++i)
    BOOST_CHECK_EQUAL(original[i], orig[i]);
}

BOOST_AUTO_TEST_SUITE_END()


} //namespace tests
} //namespace bwtc

