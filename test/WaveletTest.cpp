/**
 * @file WaveletTest.cpp
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
 * Unit tests for wavelet tree.
 */

#define BOOST_TEST_MODULE 
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include <utility>
#include <vector>

#include "../Utils.hpp"
#include "../WaveletTree.hpp"

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

BOOST_AUTO_TEST_SUITE(TreeConstructionTests)

template <typename BitVector>
int checkHuffmanShape(TreeNode<BitVector> *node, const char *answ, int *depths,
                      int curr, int depth)
{
  if(node->m_left == 0 && node->m_right == 0) {
    AlphabeticNode<BitVector>* n = (AlphabeticNode<BitVector>*) node;
    BOOST_CHECK_EQUAL(n->m_symbol, answ[curr]);
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

void checkEqual(const std::vector<bool>& vec, bool *answ) {
  for(size_t i = 0; i < vec.size(); ++i) {
    BOOST_CHECK_EQUAL(vec[i], answ[i]);
  }
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

  std::vector<bool> *codes[256];
  WaveletTree<std::vector<bool> >::collectCodes(codes, root);
  bool aCode[] = {false};
  bool bCode[] = {true, false};
  bool cCode[] = {true, true};
  checkEqual(*codes['a'], aCode);
  checkEqual(*codes['b'], bCode);
  checkEqual(*codes['c'], cCode);
  delete codes['a'];
  delete codes['b'];
  delete codes['c'];
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

  std::vector<bool> *codes[256];
  WaveletTree<std::vector<bool> >::collectCodes(codes, root);
  bool dCode[] = {false};
  bool bCode[] = {true, false, false};
  bool cCode[] = {true, false, true};
  bool aCode[] = {true, true};
  checkEqual(*codes['d'], dCode);
  checkEqual(*codes['b'], bCode);
  checkEqual(*codes['c'], cCode);
  checkEqual(*codes['a'], aCode);
  delete codes['a'];
  delete codes['b'];
  delete codes['c'];
  delete codes['d'];
}

BOOST_AUTO_TEST_CASE(HuffmanShape3) {
  uint64 freqs[256] = {0};
  const char *str = "baaabaaabcb";
  utils::calculateRunFrequencies(freqs, (const byte*) str, 11);
  TreeNode<std::vector<bool> > *root =
      WaveletTree<std::vector<bool> >::createHuffmanShape(freqs);
  const char *answers = "bac";
  int depths[] = {1, 2, 2};
  checkHuffmanShape(root, answers, depths, 0, 0);

  std::vector<bool> *codes[256];
  WaveletTree<std::vector<bool> >::collectCodes(codes, root);
  bool bCode[] = {false};
  bool aCode[] = {true, false};
  bool cCode[] = {true, true};
  checkEqual(*codes['a'], aCode);
  checkEqual(*codes['b'], bCode);
  checkEqual(*codes['c'], cCode);
  delete codes['a'];
  delete codes['b'];
  delete codes['c'];
}

BOOST_AUTO_TEST_SUITE_END()

} //namespace tests
} //namespace long_sequences

