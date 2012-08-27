/**
 * @file WaveletTest.cpp
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
 * Unit tests for wavelet tree.
 */

#define BOOST_TEST_MODULE 
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include <cstring>
#include <iterator>
#include <utility>
#include <vector>
#include <cstdlib>
#include <ctime>

#include "../Utils.hpp"
#include "../WaveletTree.hpp"
#include "TestStreams.hpp"

namespace bwtc {
int verbosity = 0;

namespace tests {

BOOST_AUTO_TEST_SUITE(HeapTests)

BOOST_AUTO_TEST_CASE(HeapTest1) {
  MinimumHeap<int> heap;
  heap.insert(4, 99);
  heap.insert(18, 3);
  heap.insert(16, 77);
  BOOST_CHECK_EQUAL(heap.deleteMin().first, 18);
  BOOST_CHECK_EQUAL(heap.deleteMin().first, 16);
  BOOST_CHECK_EQUAL(heap.deleteMin().first, 4);
}

BOOST_AUTO_TEST_CASE(HeapTest2) {
  MinimumHeap<int> heap;
  heap.insert(4, 4);
  heap.insert(20, 5);
  heap.insert(16, 16);
  heap.insert(20, 20);
  heap.insert(4, 22);
  heap.insert(16, 17);
  BOOST_CHECK_EQUAL(heap.deleteMin().first, 4);
  BOOST_CHECK_EQUAL(heap.deleteMin().first, 20);
  BOOST_CHECK_EQUAL(heap.deleteMin().first, 16);
  BOOST_CHECK_EQUAL(heap.deleteMin().first, 16);
  BOOST_CHECK_EQUAL(heap.deleteMin().first, 20);
  BOOST_CHECK_EQUAL(heap.deleteMin().first, 4);
}


BOOST_AUTO_TEST_SUITE_END()

template <typename T, typename U>
void checkEqual(const T& vec, const U& answ) {
  for(size_t i = 0; i < vec.size(); ++i) {
    BOOST_CHECK_EQUAL(vec[i], answ[i]);
  }
}

void genData(std::vector<byte>& data, int i, size_t len) {
  for(size_t j = 0; j < len; ++j) {
    int r = rand();
    // Make runs of same character more probable:
    if(j > 0 && (r & ((2 << i) - 1)) >= 1)
      data.push_back(data.back());
    else 
        data.push_back(r & 0xff);
  }
}

BOOST_AUTO_TEST_SUITE(TreeConstructionTests)

template <typename BitVector>
int checkHuffmanShape(TreeNode<BitVector> *node, const char *answ, int *depths,
                      int curr, int depth)
{
  if(node->m_left == 0 && node->m_right == 0) {
    BOOST_CHECK_EQUAL(node->m_symbol, answ[curr]);
    BOOST_CHECK_EQUAL(depth, depths[curr]);
    return curr+1;
  }
  if(node->m_left) {
    curr = checkHuffmanShape(node->m_left, answ, depths, curr, depth+1);
  }
  if(node->m_right) {
    curr = checkHuffmanShape(node->m_right, answ, depths, curr, depth+1);
  }
  return curr;
}

BOOST_AUTO_TEST_CASE(HuffmanShape1) {
  uint64 freqs[256] = {0};
  freqs['a'] = 4;
  freqs['b'] = 2;
  freqs['c'] = 1;
  TreeNode<std::vector<bool> > *root =
      WaveletTree<std::vector<bool> >::createHuffmanShape(freqs);
  const char *answers = "abc";
  int depths[] = {1,2,2};
  checkHuffmanShape(root, answers, depths, 0, 0);

  std::vector<bool> codes[256];
  WaveletTree<std::vector<bool> >::collectCodes(codes, root);
  bool aCode[] = {false};
  bool bCode[] = {true, false};
  bool cCode[] = {true, true};
  checkEqual(codes['a'], aCode);
  checkEqual(codes['b'], bCode);
  checkEqual(codes['c'], cCode);
}

BOOST_AUTO_TEST_CASE(HuffmanShape2) {
  uint64 freqs[256] = {0};
  freqs['c'] = 4;
  freqs['b'] = 5;
  freqs['a'] = 6;
  freqs['d'] = 20;
  TreeNode<std::vector<bool> > *root =
      WaveletTree<std::vector<bool> >::createHuffmanShape(freqs);
  const char *answers = "dbca";
  int depths[] = {1, 3, 3, 2};  
  checkHuffmanShape(root, answers, depths, 0, 0);

  std::vector<bool> codes[256];
  WaveletTree<std::vector<bool> >::collectCodes(codes, root);
  bool dCode[] = {false};
  bool bCode[] = {true, false, false};
  bool cCode[] = {true, false, true};
  bool aCode[] = {true, true};
  checkEqual(codes['d'], dCode);
  checkEqual(codes['b'], bCode);
  checkEqual(codes['c'], cCode);
  checkEqual(codes['a'], aCode);
}

BOOST_AUTO_TEST_CASE(HuffmanShape3) {
  uint64 freqs[256] = {0};
  const char *str = "baaabaaabcb";
  utils::calculateRunFrequencies(freqs, (const byte*) str, strlen(str));
  TreeNode<std::vector<bool> > *root =
      WaveletTree<std::vector<bool> >::createHuffmanShape(freqs);
  const char *answers = "bac";
  int depths[] = {1, 2, 2};
  checkHuffmanShape(root, answers, depths, 0, 0);

  std::vector<bool> codes[256];
  WaveletTree<std::vector<bool> >::collectCodes(codes, root);
  bool bCode[] = {false};
  bool aCode[] = {true, false};
  bool cCode[] = {true, true};
  checkEqual(codes['a'], aCode);
  checkEqual(codes['b'], bCode);
  checkEqual(codes['c'], cCode);
}

BOOST_AUTO_TEST_CASE(HuffmanShape4) {
  uint64 freqs[256] = {0};
  const char *str = "aaaa";
  utils::calculateRunFrequencies(freqs, (const byte*) str, strlen(str));
  TreeNode<std::vector<bool> > *root =
      WaveletTree<std::vector<bool> >::createHuffmanShape(freqs);
  const char *answers = "a";
  int depths[] = {1};
  checkHuffmanShape(root, answers, depths, 0, 0);

  std::vector<bool> codes[256];
  WaveletTree<std::vector<bool> >::collectCodes(codes, root);
  bool aCode[] = {false};
  checkEqual(codes['a'], aCode);
}

BOOST_AUTO_TEST_CASE(WholeConstruction1) {
  const char *str = "aaabbaaacbcb";
  WaveletTree<std::vector<bool> > tree((const byte*) str, strlen(str));
  std::vector<byte> msg;
  tree.message(std::back_inserter(msg));
  BOOST_CHECK_EQUAL(msg.size(), strlen(str));
  checkEqual(msg, (const byte*) str);
}

BOOST_AUTO_TEST_CASE(WholeConstruction2) {
  const char *str = "abbbabaagggffllslwerkfdskofdsksasdadsas"
      "dfgdfsmldsgklmesgfklmfeeeeeeeeeg";
  WaveletTree<std::vector<bool> > tree((const byte*) str, strlen(str));
  std::vector<byte> msg;
  tree.message(std::back_inserter(msg));
  BOOST_CHECK_EQUAL(msg.size(), strlen(str));
  checkEqual(msg, (const byte*) str);
}

BOOST_AUTO_TEST_CASE(WholeConstruction3) {
  const char *str = "aaaaaaaaaaaaaac";
  WaveletTree<std::vector<bool> > tree((const byte*) str, strlen(str));
  std::vector<byte> msg;
  tree.message(std::back_inserter(msg));
  BOOST_CHECK_EQUAL(msg.size(), strlen(str));
  checkEqual(msg, (const byte*) str);
}

BOOST_AUTO_TEST_CASE(WholeConstruction4) {
  const char *str = "aaaaaa"; 
  WaveletTree<std::vector<bool> > tree((const byte*) str, strlen(str));
  std::vector<byte> msg;
  tree.message(std::back_inserter(msg));
  BOOST_CHECK_EQUAL(msg.size(), strlen(str));
  checkEqual(msg, (const byte*) str);
}

BOOST_AUTO_TEST_CASE(WholeConstruction5) {
  const char *str = "abcdefghijklmnababcabcdabcdeabcd"
      "efacbcdefgabcdefghabcdefghiabcdefghij";
  WaveletTree<std::vector<bool> > tree((const byte*) str, strlen(str));
  std::vector<byte> msg;
  tree.message(std::back_inserter(msg));
  BOOST_CHECK_EQUAL(msg.size(), strlen(str));
  checkEqual(msg, (const byte*) str);
}

BOOST_AUTO_TEST_CASE(WholeConstruction6) {
  const char *str = "abaabaaabaaaabaaaaabaaaaaabaaaaaaaabaaaaaaaaaaaa";
  WaveletTree<std::vector<bool> > tree((const byte*) str, strlen(str));
  std::vector<byte> msg;
  tree.message(std::back_inserter(msg));
  BOOST_CHECK_EQUAL(msg.size(), strlen(str));
  checkEqual(msg, (const byte*) str);
}

BOOST_AUTO_TEST_CASE(WholeConstruction7) {
  char *str = new char[10000];
  for(size_t i = 0; i < 3330; ++i) {
    str[3*i] = 'a';
    str[3*i + 1] = 'b';
    str[3*i + 2] = 'c';
  }
  for(size_t i = 0; i < 4; ++i) {
    str[9990 + i] = 'a';
    str[9990 + i + 4] = 'b';
  }
  str[9998] = 'c';
  str[9999] = 'c';
  WaveletTree<std::vector<bool> > tree((const byte*) str, strlen(str));
  std::vector<byte> msg;
  tree.message(std::back_inserter(msg));
  BOOST_CHECK_EQUAL(msg.size(), strlen(str));
  checkEqual(msg, (const byte*) str);
  delete [] str;
}

BOOST_AUTO_TEST_CASE(WholeConstruction8) {
  srand(time(0));
  for(size_t i = 0; i < 5; ++i) {
    std::vector<byte> data;
    genData(data, i, 10000);
    WaveletTree<std::vector<bool> > tree(&data[0], data.size());
    std::vector<byte> msg;
    tree.message(std::back_inserter(msg));
    BOOST_CHECK_EQUAL(msg.size(), data.size());
    checkEqual(msg, data);
  }
}

BOOST_AUTO_TEST_SUITE_END()


class ByteVecWrapper {
 public:
  ByteVecWrapper() : m_bitsUsed(0), m_currByte(0) {
    m_data.push_back(0);
  }

  void push_back(bool bit) {
    int shift = 7 - m_bitsUsed;
    if(bit) {
      m_data[m_currByte] |= (1 << shift);
    } else if((m_data[m_currByte] >> shift) & 1) {
      m_data[m_currByte] &= ~(1 << shift);
    }

    ++m_bitsUsed;
    checkAndGrow();
  }

  void reset() { m_currByte = m_bitsUsed = 0; }

  byte readByte() {
    assert(m_bitsUsed == 0);
    return m_data[m_currByte++];
  }
  
  bool readBit() {
    bool ret = (m_data[m_currByte] >> (7 - m_bitsUsed)) & 1;
    ++m_bitsUsed;
    checkAndGrow();
    return ret;
  }

  std::vector<byte> m_data;

 private:
  void checkAndGrow() {
    if(m_bitsUsed == 8) {
      m_bitsUsed = 0;
      ++m_currByte;
      if(m_currByte == m_data.size()) m_data.push_back(0);
    }
  }

  uint32 m_bitsUsed;
  uint32 m_currByte;
};

BOOST_AUTO_TEST_SUITE(WaveletTreeShape)

void makeTreeShapeTest(size_t len) {
  for(size_t i = 0; i < 5; ++i) {
    std::vector<byte> data;
    genData(data, i, len);
    WaveletTree<std::vector<bool> > tree(&data[0], data.size());
    ByteVecWrapper shape;
    tree.treeShape(shape);

    shape.reset();
    
    WaveletTree<std::vector<bool> > other;
    other.readShape(shape);
    ByteVecWrapper oShape;
    other.treeShape(oShape);
    BOOST_CHECK_EQUAL(shape.m_data.size(), oShape.m_data.size());
    checkEqual(shape.m_data, oShape.m_data);
  }
}

BOOST_AUTO_TEST_CASE(TreeShape1) {
  srand(time(0));
  makeTreeShapeTest(100);
  makeTreeShapeTest(1000);
  makeTreeShapeTest(10000);
  makeTreeShapeTest(100000);
}

BOOST_AUTO_TEST_SUITE_END()

class MockCoder : public ByteVecWrapper {
 public:
  void encode(bool bit, Probability) {
    push_back(bit);
  }

  bool decode(Probability) {
    return readBit();
  }
};

class MockProbModel {
 public:
  Probability probabilityOfOne() const { return 1; }
  void update(bool) {}
  void updateState(bool) {}
  void resetModel() {}
};

BOOST_AUTO_TEST_SUITE(CodingAndDecoding)

void makeCodingAndDecodingTest(size_t len) {
  for(size_t i = 0; i < 5; ++i) {
    std::vector<byte> data;
    genData(data, i, len);
    WaveletTree<std::vector<bool> > tree(&data[0], data.size());
    MockCoder coded;
    tree.treeShape(coded);
    MockProbModel prob;
    tree.encodeTreeBF(coded, prob, prob, prob);

    size_t bitsInRoot = tree.bitsInRoot();
    coded.reset();
    
    WaveletTree<std::vector<bool> > other;
    other.readShape(coded);
    other.decodeTreeBF(bitsInRoot, coded, prob, prob, prob);

    std::vector<byte> result(data.size());
    other.message(&result[0]);
    BOOST_CHECK_EQUAL(data.size(), result.size());
    checkEqual(data, result);
  }
}

BOOST_AUTO_TEST_CASE(CodingAndDecoding1) {
  srand(time(0));
  makeCodingAndDecodingTest(100);
  makeCodingAndDecodingTest(1000);
  makeCodingAndDecodingTest(10000);
  makeCodingAndDecodingTest(100000);
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(GammaCodes)

BOOST_AUTO_TEST_CASE(Construction1) {
  std::vector<bool> one, five, seven, fifty;
  WaveletTree<std::vector<bool> >::gammaCode(one, 1);
  BOOST_CHECK_EQUAL(one.size(), 1);
  WaveletTree<std::vector<bool> >::gammaCode(five, 5);
  BOOST_CHECK_EQUAL(five.size(), 5);
  WaveletTree<std::vector<bool> >::gammaCode(seven, 7);
  BOOST_CHECK_EQUAL(seven.size(), 5);
  WaveletTree<std::vector<bool> >::gammaCode(fifty, 50);
  BOOST_CHECK_EQUAL(fifty.size(), 11);

  bool oneCode[] = {false};
  bool fiveCode[] = {true, true, false, false, true};
  bool sevenCode[] = {true, true, false, true, true};
  bool fiftyCode[] = {true, true, true, true, true,
                      false, true, false, false, true, false};
  checkEqual(one, oneCode);
  checkEqual(five, fiveCode);
  checkEqual(seven, sevenCode);
  checkEqual(fifty, fiftyCode);
}

BOOST_AUTO_TEST_SUITE_END()


} //namespace tests
} //namespace bwtc

