/**
 * @file WaveletTree.hpp
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
 * Implementation for wavelet tree.
 */

#ifndef BWTC_WAVELET_TREE_HPP_
#define BWTC_WAVELET_TREE_HPP_

#include "globaldefs.hpp"
#include "Utils.hpp"

#include <algorithm>
#include <vector>
#include <utility>

namespace bwtc {

template <typename BitVector>
struct TreeNode {
  TreeNode() : m_left(0), m_right(0) {}
  TreeNode(TreeNode<BitVector>* left, TreeNode<BitVector>* right)
      : m_left(left), m_right(right) {}

  BitVector m_bitVector;
  TreeNode<BitVector>* m_left;
  TreeNode<BitVector>* m_right;
};

// Contains character from the source alphabet
template <typename BitVector>
struct AlphabeticNode : public TreeNode<BitVector> {
  AlphabeticNode(byte symbol) : m_symbol(symbol) {}

  byte m_symbol;
};


/**This wavelet tree is used for storing the sequence
 * (<a1, n1>, <a2, n2>, ...) where a's are alphabets of the source alphabet
 * and n's are integers. Each leaf in a traditional wavelet tree is the root
 * for the wavelet tree of runs of the leaf's character. The runs (integers) are
 * coded with Elias' gamma-coding, with a minor modification.
 * For example codes for 4, 7 and 11 are the following:
 *   4 -> 11000
 *   7 -> 11011
 *  11 -> 1110011
 * The idea is to revert (log n + 1) first bits of the normal gamma code for n.
 * The original gamma-code for 11 would be 000 1 011, and in modified we revert
 * 3 + 1 first bits. This way the bit 1 means always that the coded number is
 * increasing.
 *
 * Wavelet tree is Huffman-shaped (based on the number of occurrences of the
 * runs).
 *
 * Template parameter BitVector has to implement the following member-functions:
 *   BitVector(const BitVector& other);
 *   void push_back(bool);
 *   void pop_back();
 *   void reserve(size_t size);
 *   size_t size();
 */
template <typename BitVector>
class WaveletTree {
 public:
  WaveletTree(byte *src, size_t length);
  ~WaveletTree();

  static TreeNode<BitVector> *createHuffmanShape(const uint64 *runFreqs);
  static void collectCodes(BitVector **codes, TreeNode<BitVector> *root);
  static void collectCodes(BitVector **codes, BitVector& vec, TreeNode<BitVector> *node);
  
 private:
  TreeNode<BitVector>* m_root;
  BitVector* m_codes[256];
};

template <typename BitVector>
WaveletTree<BitVector>::WaveletTree(byte *src, size_t length) {
  uint64 runFreqs[256] = {0};
  utils::calculateRunFrequencies(runFreqs, src, length);
  m_root = createHuffmanShape(runFreqs);
  collectCodes(m_codes);
}


/** This class is used on constructing Huffman-shaped wavelet tree.
 *  Essentially this is just a minimum heap which stores the weights of
 *  the nodes also. Weight for the node has to be inserted before the
 *  node itself is inserted into the heap.
 */
template <typename T>
class MinimumHeap {
 public:
  MinimumHeap() {}

  MinimumHeap(size_t initialSize) {
    m_nodes.reserve(initialSize);
  }

  void insert(T value, size_t weight) {
    size_t index = m_nodes.size();
    m_nodes.resize(m_nodes.size()+1);
    while(index > 0 && m_nodes[parent(index)].second > weight) {
      m_nodes[index] = m_nodes[parent(index)];
      index = parent(index);
    }
    m_nodes[index] = std::make_pair(value, weight);
  }

  std::pair<T, size_t> deleteMin() {
    std::pair<T, size_t> min = m_nodes[0];
    m_nodes[0] = m_nodes.back();
    m_nodes.resize(m_nodes.size()-1);
    heapify(0);
    return min;
  }

  void heapify(size_t index) {
    if(m_nodes.size() > right(index)) {
      size_t wl = m_nodes[left(index)].second;
      size_t wr = m_nodes[right(index)].second;
      size_t smaller = (wl < wr)?left(index):right(index);
      if(m_nodes[smaller].second < m_nodes[index].second) {
        std::swap(m_nodes[index], m_nodes[smaller]);
        heapify(smaller);
      } 
    } else if(m_nodes.size() == right(index) &&
              m_nodes[left(index)].second < m_nodes[index].second)
    {
      std::swap(m_nodes[index], m_nodes[left(index)]);
    }
  }

  bool empty() const { return m_nodes.size() == 0; }
  size_t size() const { return m_nodes.size(); }
  
  inline static size_t parent(size_t n) { return (n-1)/2; }
  inline static size_t left(size_t n) { return 2*n + 1; }
  inline static size_t right(size_t n) { return 2*n + 2; }

 private:
  std::vector<std::pair<T, size_t> > m_nodes;
};


template <typename BitVector>
TreeNode<BitVector>*
WaveletTree<BitVector>::createHuffmanShape(const uint64 *runFreqs) {
  MinimumHeap<TreeNode<BitVector>* > heap(256);
  byte b = 0U;
  for(int i = 0; i < 256; ++i,++b) {
    if(runFreqs[b] != 0) {
      TreeNode<BitVector> *node = new AlphabeticNode<BitVector>(b);
      heap.insert(node, runFreqs[b]);
    }
  }
  while(heap.size() > 1) {
    std::pair<TreeNode<BitVector>*, size_t> p1, p2;
    p1 = heap.deleteMin();
    p2 = heap.deleteMin();
    TreeNode<BitVector> *node = new TreeNode<BitVector>(p2.first, p1.first);
    heap.insert(node, p1.second + p2.second);
  }
  return heap.deleteMin().first;
}

template <typename BitVector>
void WaveletTree<BitVector>::collectCodes(BitVector *codes[], TreeNode<BitVector> *root) {
  BitVector vec;
  std::fill(codes, codes+256, static_cast<BitVector*>(0));
  collectCodes(codes, vec, root);
}

template <typename BitVector>
void WaveletTree<BitVector>::collectCodes(BitVector *codes[], BitVector& vec, TreeNode<BitVector> *node) {
  if(node->m_left == 0 && node->m_right == 0) {
    AlphabeticNode<BitVector>* n = (AlphabeticNode<BitVector>*) node;
    codes[n->m_symbol] = new BitVector(vec);
  }
  if(node->m_left) {
    vec.push_back(false);
    collectCodes(codes, vec, node->m_left);
    vec.pop_back();
  }
  if(node->m_right) {
    vec.push_back(true);
    collectCodes(codes, vec, node->m_right);
    vec.pop_back();
  }
}


} //namespace bwtc
#endif
