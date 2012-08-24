/**
 * @file CompressorAndDecompressorTest.cpp
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
 * Testing that Compressor and Decompressor agree on the format of compressed
 * file.
 */

#define BOOST_TEST_MODULE 
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include <cstdlib>
#include <ctime>
#include <vector>

#include "../Compressor.hpp"
#include "../Decompressor.hpp"
#include "TestStreams.hpp"

namespace bwtc {
int verbosity = 0;

namespace tests {

void makeRandomData(std::vector<byte>& data, size_t length) {
  for(size_t i = 0; i < length; ++i) {
    byte b = rand() & 0xff;
    data.push_back(b);
  }
}

void makeRepetitiveData(std::vector<byte>& data, size_t length, int times) {
  makeRandomData(data, length);
  for(int j = 0; j < times-1; ++j) {
    for(size_t i = 0; i < length; ++i) {
        data.push_back(data[i]);
    }
  }
}

void test(size_t length, size_t reps, const char* prep, size_t mem,
          char entropyCoder, char bwtAlgo, size_t startingPoints)
{
  srand(time(0));
  std::vector<byte> orig, comp, decomp;
  TestStream *original = new TestStream(orig),
      *compressed = new TestStream(comp),
      *compr2 = new TestStream(comp),
      *decompressed = new TestStream(decomp);
  if(reps == 0) {
    makeRandomData(orig, length);
  } else {
    makeRepetitiveData(orig, length/reps, reps);
  }

  {
    Compressor compressor(original, compressed, prep, mem,
                          entropyCoder);
    compressor.initializeBwtAlgorithm(bwtAlgo, startingPoints);
    compressor.compress(1);
    
    Decompressor decompressor(compr2, decompressed);
    decompressor.decompress(1);

    BOOST_CHECK_EQUAL(orig.size(), decomp.size());
    for(size_t i = 0; i < orig.size(); ++i)
      BOOST_CHECK_EQUAL(orig[i], decomp[i]);
  }
}


BOOST_AUTO_TEST_SUITE(WithWaveletCoders)

BOOST_AUTO_TEST_CASE(SingleStartingPointSingleBlockNoPreprocessing) {
  test(100, 0, "", 1000, 'B', 'd', 1);
  test(1000, 0, "", 10000, 'B', 'd', 1);
  test(10000, 0, "", 100000, 'B', 'd', 1);
  test(100000, 0, "", 1000000, 'B', 'd', 1);

  test(100, 2, "", 1000, 'B', 's', 1);
  test(1000, 2, "", 10000, 'B', 's', 1);
  test(10000, 2, "", 100000, 'B', 's', 1);
  test(100000, 2, "", 1000000, 'B', 's', 1);

  test(100, 50, "", 1000, 'B', 'd', 1);
  test(1000, 50, "", 10000, 'B', 'd', 1);
  test(10000, 50, "", 100000, 'B', 'd', 1);
  test(100000, 50, "", 1000000, 'B', 'd', 1);
}

BOOST_AUTO_TEST_CASE(SingleStartingPointSingleBlockPreprocessing) {
  test(100, 0, "p", 1000, 'B', 'd', 1);
  test(1000, 0, "p", 10000, 'B', 'd', 1);
  test(10000, 0, "p", 100000, 'B', 'd', 1);
  test(100000, 0, "p", 1000000, 'B', 'd', 1);

  test(100, 2, "pp", 1000, 'B', 's', 1);
  test(1000, 2, "pp", 10000, 'B', 's', 1);
  test(10000, 2, "pp", 100000, 'B', 's', 1);
  test(100000, 2, "pp", 1000000, 'B', 's', 1);

  test(100, 50, "pppp", 1000, 'B', 'd', 1);
  test(1000, 50, "ppp", 10000, 'B', 'd', 1);
  test(10000, 50, "ppp", 100000, 'B', 'd', 1);
  test(100000, 50, "ppppp", 1000000, 'B', 'd', 1);
}

BOOST_AUTO_TEST_CASE(MultipleStartingPoints) {
  for(size_t i = 1; i <= 30; ++i) 
    test(10000, 0, "", 100000, 'B', 'd', i);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(WithHuffmanCoders)

BOOST_AUTO_TEST_CASE(SingleStartingPointSingleBlockNoPreprocessing) {
  test(100, 0, "", 1000, 'H', 'd', 1);
  test(1000, 0, "", 10000, 'H', 'd', 1);
  test(10000, 0, "", 100000, 'H', 'd', 1);
  test(100000, 0, "", 1000000, 'H', 'd', 1);

  test(100, 2, "", 1000, 'H', 's', 1);
  test(1000, 2, "", 10000, 'H', 's', 1);
  test(10000, 2, "", 100000, 'H', 's', 1);
  test(100000, 2, "", 1000000, 'H', 's', 1);

  test(100, 50, "", 1000, 'H', 'd', 1);
  test(1000, 50, "", 10000, 'H', 'd', 1);
  test(10000, 50, "", 100000, 'H', 'd', 1);
  test(100000, 50, "", 1000000, 'H', 'd', 1);
}

BOOST_AUTO_TEST_CASE(SingleStartingPointSingleBlockPreprocessing) {
  test(100, 0, "p", 1000, 'H', 'd', 1);
  test(1000, 0, "p", 10000, 'H', 'd', 1);
  test(10000, 0, "p", 100000, 'H', 'd', 1);
  test(100000, 0, "p", 1000000, 'H', 'd', 1);

  test(100, 2, "pp", 1000, 'H', 's', 1);
  test(1000, 2, "pp", 10000, 'H', 's', 1);
  test(10000, 2, "pp", 100000, 'H', 's', 1);
  test(100000, 2, "pp", 1000000, 'H', 's', 1);

  test(100, 50, "pppp", 1000, 'H', 'd', 1);
  test(1000, 50, "ppp", 10000, 'H', 'd', 1);
  test(10000, 50, "ppp", 100000, 'H', 'd', 1);
  test(100000, 50, "ppppp", 1000000, 'H', 'd', 1);
}

BOOST_AUTO_TEST_CASE(MultipleStartingPoints) {
  for(size_t i = 1; i <= 30; ++i) 
    test(10000, 0, "", 100000, 'H', 'd', i);
}

BOOST_AUTO_TEST_SUITE_END()


} //namespace tests
} //namespace bwtc

