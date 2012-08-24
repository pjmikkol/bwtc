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

#include <vector>

#include "../globaldefs.hpp"
#include "../Utils.hpp"

#include "TestStreams.hpp"

namespace bwtc {
int verbosity = 0;

namespace tests {

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

BOOST_AUTO_TEST_CASE(LogFlrU) {
  BOOST_CHECK_EQUAL(7, logFloor(0xffu));
  BOOST_CHECK_EQUAL(26, logFloor(0x61adf3fu));
  BOOST_CHECK_EQUAL(30, logFloor(0x561adf3fu));
  BOOST_CHECK_EQUAL(31, logFloor(0xf61adf3fu));
}

BOOST_AUTO_TEST_CASE(LogFlrUL) {
  BOOST_CHECK_EQUAL(7, logFloor(0xfful));
  BOOST_CHECK_EQUAL(26, logFloor(0x61adf3ful));
  BOOST_CHECK_EQUAL(30, logFloor(0x561adf3ful));
  BOOST_CHECK_EQUAL(31, logFloor(0xf61adf3ful));
}

BOOST_AUTO_TEST_CASE(RunFrequencyCounts1) {
  const char * str = "aabacca";
  uint64 freqs[256] = {0};
  calculateRunFrequencies(freqs, (const byte*)str, strlen(str));
  BOOST_CHECK_EQUAL(freqs['a'], 3);
  BOOST_CHECK_EQUAL(freqs['c'], 1);
  BOOST_CHECK_EQUAL(freqs['b'], 1);
  for(int i = 0; i < 256; ++i) {
    if(i == 'a' || i == 'c' || i == 'b') continue;
    BOOST_CHECK_EQUAL(freqs[i],0);
  }
}

BOOST_AUTO_TEST_CASE(RunFrequencyCounts2) {
  const char * str = "wwwaaaawawawab";
  uint64 freqs[256] = {0};
  calculateRunFrequencies(freqs, (const byte*)str, strlen(str));
  BOOST_CHECK_EQUAL(freqs['w'], 4);
  BOOST_CHECK_EQUAL(freqs['a'], 4);
  BOOST_CHECK_EQUAL(freqs['b'], 1);
  for(int i = 0; i < 256; ++i) {
    if(i == 'a' || i == 'w' || i == 'b') continue;
    BOOST_CHECK_EQUAL(freqs[i],0);
  }
}

BOOST_AUTO_TEST_CASE(RunFrequencyCounts3) {
  const char * str = "baaabaaabcb";
  uint64 freqs[256] = {0};
  calculateRunFrequencies(freqs, (const byte*)str, strlen(str));
  BOOST_CHECK_EQUAL(freqs['a'], 2);
  BOOST_CHECK_EQUAL(freqs['b'], 4);
  BOOST_CHECK_EQUAL(freqs['c'], 1);
  for(int i = 0; i < 256; ++i) {
    if(i == 'a' || i == 'c' || i == 'b') continue;
    BOOST_CHECK_EQUAL(freqs[i],0);
  }
}

BOOST_AUTO_TEST_CASE(RunFrequencyCounts4) {
  const char * str = "baaabaaabcbb";
  uint64 freqs[256] = {0};
  calculateRunFrequencies(freqs, (const byte*)str, strlen(str));
  BOOST_CHECK_EQUAL(freqs['a'], 2);
  BOOST_CHECK_EQUAL(freqs['b'], 4);
  BOOST_CHECK_EQUAL(freqs['c'], 1);
  for(int i = 0; i < 256; ++i) {
    if(i == 'a' || i == 'c' || i == 'b') continue;
    BOOST_CHECK_EQUAL(freqs[i],0);
  }
}


#define INTEGER_PACKING_TEST(par) \
  std::vector<byte> vec;\
  TestStream stream(vec); \
  int bytesWritten; \
  uint64 integer = (par);\
  uint64 packed = utils::packInteger(integer, &bytesWritten);\
  for(int i = 0; i < bytesWritten; ++i) {\
    stream.writeByte(packed & 0xff);\
    packed >>= 8;\
  }\
  stream.reset();\
  size_t bytesRead;\
  size_t result = utils::readPackedInteger(stream, bytesRead);\
  BOOST_CHECK_EQUAL(result, integer);\
  BOOST_CHECK_EQUAL(bytesRead, bytesWritten)

BOOST_AUTO_TEST_CASE(PackedInteger) {
  uint64 numbers[] = {0, 1};//, 12, 126, 127, 128, 250, 1234, 5422, 124312, 4311235};
  for(size_t i = 0; i < sizeof(numbers)/sizeof(uint64); ++i) {
    INTEGER_PACKING_TEST(numbers[i]);
  }
}


BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(HuffmanLengths)

BOOST_AUTO_TEST_CASE(HuffmanLengths1) {
  uint64 freqs[256] = {0};
  freqs['a'] = 5;
  freqs['b'] = 1;
  freqs['c'] = 1;
  std::vector<std::pair<uint64, uint32> > t;
  calculateHuffmanLengths(t, freqs);
  BOOST_CHECK_EQUAL(t[2].first, 1);
  BOOST_CHECK_EQUAL(t[2].second, 'a');
  BOOST_CHECK_EQUAL(t[0].first, 2);
  BOOST_CHECK_EQUAL(t[0].second, 'b');
  BOOST_CHECK_EQUAL(t[1].first, 2);
  BOOST_CHECK_EQUAL(t[1].second, 'c');
}

BOOST_AUTO_TEST_CASE(HuffmanLengths2) {
  uint64 freqs[256] = {0};
  freqs['a'] = 1;
  freqs['b'] = 3;
  freqs['c'] = 4;
  freqs['d'] = 6;
  freqs['e'] = 8;
  freqs['f'] = 4;
  freqs['g'] = 1;
  std::vector<std::pair<uint64, uint32> > t;
  calculateHuffmanLengths(t, freqs);
  BOOST_CHECK_EQUAL(t[0].first, 4);
  BOOST_CHECK_EQUAL(t[1].first, 4);
  BOOST_CHECK_EQUAL(t[2].first, 3);
  BOOST_CHECK_EQUAL(t[3].first, 3);
  BOOST_CHECK_EQUAL(t[4].first, 3);
  BOOST_CHECK_EQUAL(t[5].first, 2);
  BOOST_CHECK_EQUAL(t[6].first, 2);
}

BOOST_AUTO_TEST_CASE(HuffmanLengths3) {
  uint64 freqs[256] = {0};
  freqs['a'] = 2;
  freqs['b'] = 2;
  freqs['c'] = 4;
  freqs['d'] = 8;
  freqs['e'] = 8;
  std::vector<std::pair<uint64, uint32> > t;
  calculateHuffmanLengths(t, freqs);

  BOOST_CHECK_EQUAL(t[0].first, 3);
  BOOST_CHECK_EQUAL(t[1].first, 3);
  BOOST_CHECK_EQUAL(t[2].first, 2);
  BOOST_CHECK_EQUAL(t[3].first, 2);
  BOOST_CHECK_EQUAL(t[4].first, 2);
}

BOOST_AUTO_TEST_CASE(HuffmanLengths4) {
  uint64 freqs[256] = {0};
  freqs['a'] = 2;
  freqs['b'] = 20000;
  std::vector<std::pair<uint64, uint32> > t;
  calculateHuffmanLengths(t, freqs);

  BOOST_CHECK_EQUAL(t[0].first, 1);
  BOOST_CHECK_EQUAL(t[1].first, 1);
}

BOOST_AUTO_TEST_CASE(HuffmanLengths5) {
  uint64 freqs[256] = {0};
  freqs['a'] = 1;
  std::vector<std::pair<uint64, uint32> > t;
  calculateHuffmanLengths(t, freqs);

  BOOST_CHECK_EQUAL(t[0].first, 1);
}


BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(InterpolativeCoding)


struct Input {
  Input() : curr(0) {}
  Input(const std::vector<bool>& v) : curr(0), data(v) {}
  bool readBit() { return data.at(curr++); }
  void pb(bool bit) { data.push_back(bit); }
  size_t curr; 
  std::vector<bool> data;
};


template <typename T1, typename T2>
void checkEq(const T1& t1, const T2& t2, size_t len) {
  for(size_t i = 0; i < len; ++i)  {
    BOOST_CHECK_EQUAL(t1[i],t2[i]);
  }
}

#define BINARY_CODE_TEST(N, LO, HI, ...) {   \
  std::vector<bool> bits; \
  binaryCode(N, LO, HI, bits); \
  bool t1[] = {__VA_ARGS__};          \
  size_t n = sizeof(t1)/sizeof(bool); \
  BOOST_CHECK_EQUAL(bits.size(), n); \
  checkEq(bits, t1, n); \
}

BOOST_AUTO_TEST_CASE(BinaryCode1) {
  BINARY_CODE_TEST(8, 0, 8, false, false, false, true);
  BINARY_CODE_TEST(2, 0, 8, false, true, false);
  BINARY_CODE_TEST(5, 0, 8, true, false, true);
  BINARY_CODE_TEST(0, 0, 8, false, false, false, false);
  BINARY_CODE_TEST(1, 0, 8, false, false, true);
}

BOOST_AUTO_TEST_CASE(BinaryCode2) {
  BINARY_CODE_TEST(3, 0, 14, false, false, true, true);
  BINARY_CODE_TEST(7, 0, 14, true, true, true);
  BINARY_CODE_TEST(11, 0, 14, true, false, true, false);
  BINARY_CODE_TEST(13, 0, 14, true, true, false, false);
  BINARY_CODE_TEST(5, 0, 14, false, true, false, true);
}

BOOST_AUTO_TEST_CASE(BinaryCode3) {
  BINARY_CODE_TEST(3, 0, 7, false, true, true);
  BINARY_CODE_TEST(4, 0, 7, true, false, false);
  BINARY_CODE_TEST(7, 0, 7, true, true, true);
  BINARY_CODE_TEST(14, 0, 15, true, true, true, false);
  BINARY_CODE_TEST(7, 0, 15, false, true, true, true);
}

BOOST_AUTO_TEST_CASE(BinaryCode4) {
  BINARY_CODE_TEST(10, 3, 16, true, true, true);
  BINARY_CODE_TEST(7, 1, 9, true, true, false);
  BINARY_CODE_TEST(2, 0, 6, false, true, false);
  BINARY_CODE_TEST(8, 8, 9, false);
  BINARY_CODE_TEST(12, 12, 18, false, false, false);
  BINARY_CODE_TEST(16, 13, 19, true, true);
}

#define BINARY_DECODE_TEST(N, LO, HI, ...) { \
  Input in; \
  bool t[] = {__VA_ARGS__}; \
  for(size_t i = 0; i < sizeof(t)/sizeof(bool); ++i) in.pb(t[i]); \
  BOOST_CHECK_EQUAL(binaryDecode(in, LO, HI), N); \
  BOOST_CHECK_EQUAL(in.data.size(), in.curr); \
}


BOOST_AUTO_TEST_CASE(BinaryDecode1) {
  BINARY_DECODE_TEST(8, 0, 8, false, false, false, true);
  BINARY_DECODE_TEST(2, 0, 8, false, true, false);
  BINARY_DECODE_TEST(5, 0, 8, true, false, true);
  BINARY_DECODE_TEST(0, 0, 8, false, false, false, false);
  BINARY_DECODE_TEST(1, 0, 8, false, false, true);
}

BOOST_AUTO_TEST_CASE(BinaryDecode2) {
  BINARY_DECODE_TEST(3, 0, 14, false, false, true, true);
  BINARY_DECODE_TEST(7, 0, 14, true, true, true);
  BINARY_DECODE_TEST(11, 0, 14, true, false, true, false);
  BINARY_DECODE_TEST(13, 0, 14, true, true, false, false);
  BINARY_DECODE_TEST(5, 0, 14, false, true, false, true);
}

BOOST_AUTO_TEST_CASE(BinaryDecode3) {
  BINARY_DECODE_TEST(3, 0, 7, false, true, true);
  BINARY_DECODE_TEST(4, 0, 7, true, false, false);
  BINARY_DECODE_TEST(7, 0, 7, true, true, true);
  BINARY_DECODE_TEST(14, 0, 15, true, true, true, false);
  BINARY_DECODE_TEST(7, 0, 15, false, true, true, true);
}

BOOST_AUTO_TEST_CASE(BinaryDecode4) {
  BINARY_DECODE_TEST(10, 3, 16, true, true, true);
  BINARY_DECODE_TEST(7, 1, 9, true, true, false);
  BINARY_DECODE_TEST(2, 0, 6, false, true, false);
  BINARY_DECODE_TEST(8, 8, 9, false);
  BINARY_DECODE_TEST(12, 12, 18, false, false, false);
  BINARY_DECODE_TEST(16, 13, 19, true, true);
}


#undef BINARY_CODE_TEST
#undef BINARY_DECODE_TEST
#define PACK_VECTOR(V, T, ...) { \
  T t[] = {__VA_ARGS__}; \
  for(size_t i = 0; i < sizeof(t)/sizeof(T); ++i) \
    V.push_back(t[i]); \
}

BOOST_AUTO_TEST_CASE(InterpolativeCode1) {
  std::vector<size_t> list;
  PACK_VECTOR(list, size_t, 2, 7, 8, 10, 11, 12, 16);
  std::vector<bool> code;
  binaryInterpolativeCode(list, 19, code);
  std::vector<bool> correct;
  PACK_VECTOR(correct, bool,
              true, true, true,      /* 10*/
              true, true, false,     /* 7 */
              false, true, false,    /* 2 */
              false,                 /* 8 */
              false, false, false,   /* 12*/
                                     /* 11*/
              true, true             /* 3 */
              );
  checkEq(code, correct, 15);
}

BOOST_AUTO_TEST_CASE(InterpolativeCode2) {
  std::vector<size_t> list;
  PACK_VECTOR(list, size_t, 3, 4, 6, 7);
  std::vector<bool> code;
  binaryInterpolativeCode(list, 7, code);
  std::vector<bool> correct;
  PACK_VECTOR(correct, bool,
              true, true,            /* 4 */
              true, true,            /* 3 */
              true,                  /* 6 */
                                     /* 7 */
              );
  checkEq(code, correct, 5);
}

BOOST_AUTO_TEST_CASE(InterpolativeDecode1) {
  std::vector<size_t> correct;
  PACK_VECTOR(correct, size_t, 2, 7, 8, 10, 11, 12, 16);
  std::vector<bool> code;
  PACK_VECTOR(code, bool,
              true, true, true,      /* 10*/
              true, true, false,     /* 7 */
              false, true, false,    /* 2 */
              false,                 /* 8 */
              false, false, false,   /* 12*/
                                     /* 11*/
              true, true             /* 3 */
              );
  Input in(code);
  std::vector<size_t> list;
  binaryInterpolativeDecode(list, in, 19, 7);
  checkEq(list, correct, 7);
}

BOOST_AUTO_TEST_CASE(InterpolativeDecode2) {
  std::vector<size_t> correct;
  PACK_VECTOR(correct, size_t, 3, 4, 6, 7);
  std::vector<bool> code;
  PACK_VECTOR(code, bool,
              true, true,            /* 4 */
              true, true,            /* 3 */
              true,                  /* 6 */
                                     /* 7 */
              );
  Input in(code);
  std::vector<size_t> list;
  binaryInterpolativeDecode(list, in, 7, 4);
  BOOST_CHECK_EQUAL(in.curr, 5);
  checkEq(list, correct, 4);
}


BOOST_AUTO_TEST_SUITE_END()

} //namespace tests
} //namespace long_sequences

