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
 * Implementation for wavelet tree. This variant doesn't have
 * access(i)-operation because runs of the same character are compressed
 * and the length of each run is stored into the structure of the tree.
 * On the other hand question 'which is the i-th run' can be answered
 * efficiently provided that the BitVector used has efficient implementation
 * of rank operation.
 */

#ifndef BWTC_WAVELET_TREE_HPP_
#define BWTC_WAVELET_TREE_HPP_

#include "globaldefs.hpp"
#include "Utils.hpp"
#include "Profiling.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <list>
#include <map>
#include <queue>
#include <utility>
#include <vector>

//#define OPTIMIZED_INTEGER_CODE
//#define SEMI_FIXED_CODE

namespace bwtc {

template <typename BitVector>
struct TreeNode {
  TreeNode() : m_left(0), m_right(0), m_hasSymbol(false) {}
  TreeNode(uint32 symbol) : m_left(0), m_right(0), m_symbol(symbol),
                            m_hasSymbol(true) {}
  TreeNode(TreeNode<BitVector>* left, TreeNode<BitVector>* right)
      : m_left(left), m_right(right), m_hasSymbol(false) {}
  ~TreeNode() {}

  size_t rank(bool bit, size_t i) const;
  size_t totalBits() const;

  BitVector m_bitVector;
  TreeNode<BitVector>* m_left;
  TreeNode<BitVector>* m_right;
  uint32 m_symbol;
  bool m_hasSymbol;
  
};

/**Generic rank-query for tree node. Doesn't rely on special properties
 * of BitVector and hence takes linear time.
 * If ever using something suitable for fast rank-queries, there should be
 * template specializations for rank-queries.
 */
template <typename BitVector>
size_t TreeNode<BitVector>::rank(bool bit, size_t i) const {
  if(m_bitVector.size() == 0) return 0;
  i = std::min(m_bitVector.size(), i);
  size_t sum = 0;
  for(size_t j = 0; j < i; ++j) if(!(m_bitVector[j] ^ bit)) ++sum;
  return sum;
}

template <typename BitVector>
size_t TreeNode<BitVector>::totalBits() const {
  size_t bits = m_bitVector.size();
  if(m_left) bits += m_left->totalBits();
  if(m_right) bits += m_right->totalBits();
  return bits;
}

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
 * runs). In this implementation the wavelet tree has always at least two
 * nodes (root and a leaf).
 *
 * Template parameter BitVector has to implement the following member-functions:
 *   BitVector();
 *   BitVector& operator=(const BitVector&);
 *   void push_back(bool);
 *   void pop_back();
 *   void reserve(size_t size);
 *   size_t size();
 *   bool operator[](size_t index);
 */
template <typename BitVector>
class WaveletTree {
 public:
  WaveletTree();
  WaveletTree(const byte *src, size_t length);
  ~WaveletTree();
  void pushRun(byte symbol, size_t runLength);
  TreeNode<BitVector> *pushBits(const BitVector& bits);

  /**Writes encoding of the shape of tree into given vector. The shape of
   * the tree is in form presented in "Housekeeping for Prefix Codes" by
   * Turpin and Moffat. */
  void treeShape(BitVector& vector) const;

  size_t bitsInRoot() const { return m_root->m_bitVector.size(); }
  size_t totalBits() const { return m_root->totalBits(); }

  template <typename OutputIterator>
  size_t message(OutputIterator out) const;

  /** Push the runs of the string into tree */
  void pushMessage(const byte* src, size_t length);

  /**Input is required to have readBit()-method returning bool and readByte().
   */
  template<typename Input>
  size_t readShape(Input& input);

  /** Encoder is required to have encode(bool bit, Probability prob)-method.
    * ProbabilisticModel is required to have probabilityOfOne()- and
    * update(bool bit)-methods. */
  template <typename Encoder, typename ProbabilisticModel>
  void encodeTree(Encoder& enc, ProbabilisticModel& pm) const;

  template <typename Encoder, typename ProbabilisticModel, typename GammaModel,
            typename GapModel>
#ifdef ENTROPY_PROFILER
  void encodeTreeBF(Encoder& enc, ProbabilisticModel& pm, GammaModel& gm,
                    GapModel& gapm);
#else
  void encodeTreeBF(Encoder& enc, ProbabilisticModel& pm, GammaModel& gm,
                    GapModel& gapm) const;
#endif

  /** Decoder is required to have decode(Probability prob)-method.
    * ProbabilisticModel is required to have probabilityOfOne()- and
    * update(bool bit)-methods. */
  template <typename Decoder, typename ProbabilisticModel>
  void decodeTree(size_t rootSize, Decoder& enc, ProbabilisticModel& pm) const;

  template <typename Decoder, typename ProbabilisticModel, typename GammaModel,
            typename GapModel>
  void decodeTreeBF(size_t rootSize, Decoder& enc, ProbabilisticModel& pm,
                    GammaModel& gm, GapModel& gapm) const;

  const BitVector& code(byte symbol) { return m_codes[symbol]; }
  
  static TreeNode<BitVector> *createHuffmanShape(const uint64 *runFreqs);

  template <typename BVectors>
  static void collectCodes(BVectors& codes, TreeNode<BitVector> *root);

  template <typename BVectors>
  static void collectCodes(BVectors& codes, BitVector& vec,
                           TreeNode<BitVector> *node);

  static void pushBits(TreeNode<BitVector> *node, const BitVector& bits);
  static void pushBits(TreeNode<BitVector> *node, const BitVector& bits, uint32 symbol);
  static void gammaCode(BitVector& bits, size_t integer);

  static void fixedIntegerCode(BitVector& bits, uint32 integer, uint32 w);
  static size_t fixedIntegerCodeTranslation(size_t bits, uint32 w);

  /**Assigns prefix codes based on the given lengths. Actual codes are stored
   * into m_codes.
   *
   * @param lengths Lengths of the codewords to be assigned given as
   *                <length, symbol>-pairs. This array WILL be modified.
   */
  void assignPrefixCodes(std::vector<std::pair<uint64, uint32> >& lengths);

#ifdef ENTROPY_PROFILER
  uint32 m_bytesForCharacters;
  uint32 m_bytesForRuns;
#endif

  
 private:
  TreeNode<BitVector>* m_root;
  BitVector m_codes[256];

  // Parameter of fixed integer codes
  static const size_t m_W = 6;

#ifdef OPTIMIZED_INTEGER_CODE
  std::map<uint32, BitVector> m_integerCodes;
  // used for tracking the code in decoder
  TreeNode<BitVector>* m_integerCodeTree;
#endif
  
  template <typename Decoder, typename ProbabilisticModel>
  void decodeTree(TreeNode<BitVector> *node, size_t nodeSize, Decoder& dec,
                  ProbabilisticModel& pm) const;

  template <typename Decoder, typename ProbabilisticModel>
  void decodeTreeGammaInc(TreeNode<BitVector> *node, size_t nodeSize, size_t bitsInCode,
                          Decoder& dec, ProbabilisticModel& pm) const;

  template <typename Decoder, typename ProbabilisticModel>
  void decodeTreeGammaDec(TreeNode<BitVector> *node, size_t nodeSize, size_t bitsInCode,
                          Decoder& dec, ProbabilisticModel& pm) const;

  template <typename Encoder, typename ProbabilisticModel>
  void encodeTree(TreeNode<BitVector>* node, Encoder& enc,
                  ProbabilisticModel& pm) const;

  static void destroy(TreeNode<BitVector>* node);

  /**Recursive function which is used in assigning codes and constructing the
   * tree based on the lengths of codes.
   *
   * @param lengths Lengths of codes and corresponding symbols.
   * @param node Node where the traversal is at the moment.
   * @param elem Element (index to lengths) which is handled at the moment.
   * @param bits How many bits are appended into the current element.
   * @return Next element (index) to handle.
   */
  static size_t assignPrefixCodes(std::vector<std::pair<uint64, uint32> >& lengths,
                                  TreeNode<BitVector>* node, size_t elem, size_t bits);

  /**Finds good parameters for semi-fixed coding.
   *
   * @param integerFrequencies List of integers and their frequencies
   *        (first member) is the frequency.
   * @param totalFrequencies Sum of the first members in the integerFrequencies-
   *        vector.
   * @return The pair returned contains value of W used in the computation of
   *         semi-fixed codes and the depth of the node representing the
   *         semi-fixed codes in the code-tree. 
   */
  static uint32 findParametersForSemiFixedCodes(
      std::vector<std::pair<uint64, uint32> >& integerFrequencies,
      size_t totalFrequencies);
};

template <typename BitVector>
WaveletTree<BitVector>::WaveletTree()
#ifdef ENTROPY_PROFILER
    : m_bytesForCharacters(0), m_bytesForRuns(0)
#endif
{
#ifdef OPTIMIZED_INTEGER_CODE
  m_integerCodeTree = 0;
#endif
  m_root = new TreeNode<BitVector>();
}

#ifdef OPTIMIZED_INTEGER_CODE

template <typename BitVector>
WaveletTree<BitVector>::WaveletTree(const byte *src, size_t length)
#ifdef ENTROPY_PROFILER
    : m_bytesForCharacters(0), m_bytesForRuns(0)
#endif
{
  PROFILE("WaveletTree::WaveletTree");
  uint64 runFreqs[256] = {0};

  // Frequencies of run lengths are collected into map. They are indexed
  // by the length of run.
  std::map<uint32, uint32> runDistribution;
  size_t totalRuns = utils::calculateRunsAndCharacters(
      runFreqs, src, length, runDistribution);

  // Calculate codes for the byte-alphabet (top part of Wavelet-tree)
  std::vector<std::pair<uint64, uint32> > codeLengths;
  utils::calculateHuffmanLengths(codeLengths, runFreqs);

  assignPrefixCodes(codeLengths);

  {
    assert(runDistribution.size() > 0);
    std::vector<std::pair<uint64, uint32> > integerCodeLengths;
    std::vector<uint64> freqs;
    std::vector<uint32> integers;
    for(std::map<uint32, uint32>::const_iterator it = runDistribution.begin();
        it != runDistribution.end(); ++it) {
#ifdef SEMI_FIXED_CODE
      integerCodeLengths.push_back(std::make_pair(it->second, it->first));
#else
      integers.push_back(it->first);
      freqs.push_back(it->second);
#endif
    }

#ifdef SEMI_FIXED_CODE
    //m_W = findParametersForSemiFixedCodes(integerCodeLengths, totalRuns);
    findParametersForSemiFixedCodes(integerCodeLengths, totalRuns);

    freqs.resize(integerCodeLengths.size());
    utils::calculateCodeLengths(integerCodeLengths, &freqs[0], true);
    std::sort(integerCodeLengths.begin(), integerCodeLengths.end());
#else
    /* Huffman-codes for integers */
    utils::calculateHuffmanLengths(integerCodeLengths, &freqs[0], integers);
    std::sort(integerCodeLengths.begin(), integerCodeLengths.end());
    
    /* Hu-Tucker-codes for integers
    utils::calculateHuTuckerLengths(integerCodeLengths, &freqs[0], integers);
    */
#endif
    m_integerCodeTree = new TreeNode<BitVector>();
    assignPrefixCodes(integerCodeLengths, m_integerCodeTree, 0, 0);
    collectCodes(m_integerCodes, m_integerCodeTree);
  }
  
  collectCodes(m_codes, m_root);

  pushMessage(src, length);
}

#else //ifndef OPTIMIZED_INTEGER_CODE

template <typename BitVector>
WaveletTree<BitVector>::WaveletTree(const byte *src, size_t length)
#ifdef ENTROPY_PROFILER
    : m_bytesForCharacters(0), m_bytesForRuns(0)
#endif
{
  PROFILE("WaveletTree::WaveletTree");
  uint64 runFreqs[256] = {0};

  utils::calculateRunFrequencies(runFreqs, src, length);

  // Calculate codes for the byte-alphabet (top part of Wavelet-tree)
  std::vector<std::pair<uint64, uint32> > codeLengths;
  utils::calculateHuffmanLengths(codeLengths, runFreqs);

  assignPrefixCodes(codeLengths);
  
  collectCodes(m_codes, m_root);

  pushMessage(src, length);
}
#endif

template <typename BitVector>
WaveletTree<BitVector>::~WaveletTree() {
  destroy(m_root);
#ifdef OPTIMIZED_INTEGER_CODE
  destroy(m_integerCodeTree);
#endif  
}

template <typename BitVector>
void WaveletTree<BitVector>::destroy(TreeNode<BitVector>* node) {
  if(node->m_left) destroy(node->m_left);
  if(node->m_right) destroy(node->m_right);
  delete node;
}

template <typename BitVector> template <typename Input>
size_t WaveletTree<BitVector>::readShape(Input& input) {
  assert(m_root->m_left == 0 && m_root->m_right == 0);
  size_t maxSym = input.readByte();
  size_t symbols = input.readByte();
  if(symbols == 0) symbols = 256;

  size_t bitsRead = 0;
  size_t maxLen = utils::readPackedIntegerRev(input, bitsRead);
  bitsRead = 8*bitsRead + 16;

  std::vector<byte> alphabet;
  bitsRead += utils::binaryInterpolativeDecode(alphabet, input,
                                               maxSym, symbols);

  std::vector<std::pair<uint64, uint32> > codeLengths;
  for(size_t i = 0; i < symbols; ++i) {
    size_t n = utils::unaryDecode(input);
    bitsRead += n;
    size_t len = maxLen - n + 1;
    codeLengths.push_back(std::make_pair(len, alphabet[i]));
  }

  assignPrefixCodes(codeLengths);
  collectCodes(m_codes, m_root);
  
#ifdef OPTIMIZED_INTEGER_CODE

  {
    size_t bytesRead = 0;
    size_t longestRun = utils::readPackedIntegerRev(input, bytesRead);
    bitsRead += 8*bytesRead;

    symbols = utils::readPackedIntegerRev(input, bytesRead);
    bitsRead += 8*bytesRead;


    std::vector<uint32> integers;
    maxLen = utils::readPackedIntegerRev(input, bytesRead);
    bitsRead += 8*bytesRead;
#ifndef SEMI_FIXED_CODE
    bitsRead += utils::binaryInterpolativeDecode(integers, input, 1,
                                                 longestRun, symbols);
#else
    bitsRead += utils::binaryInterpolativeDecode(integers, input, 0,
                                                 longestRun, symbols);
#endif

    std::vector<std::pair<uint64, uint32> > integerCodeLengths;
    for(size_t i = 0; i < symbols; ++i) {
      size_t n = utils::unaryDecode(input);
      bitsRead += n;
      assert(n <= maxLen);
      size_t len = maxLen + 1 - n;
      integerCodeLengths.push_back(std::make_pair(len, integers[i]));
    }

    /* Would sort if used Huffman codes*/
    std::sort(integerCodeLengths.begin(), integerCodeLengths.end());
    m_integerCodeTree = new TreeNode<BitVector>();

    assignPrefixCodes(integerCodeLengths, m_integerCodeTree, 0, 0);
    collectCodes(m_integerCodes, m_integerCodeTree);

    assert(!m_integerCodeTree->m_hasSymbol);
  }
#endif
  return bitsRead;
}

/**This code resembles gamma code, but the tree which represents it, is
 * flatter. Let W >= 0 (the bigger W, the flatter the tree is).
 * Let x be positive integer. The code C(x) of x has two parts:
 *
 * The first part tells how long the second part is. It has 
 *      B = floor( log[x - 1 + 2^W] ) - W
 * one-bits followed by single zero-bit. The logarithm has base of two.
 *
 * The second part has W + B least significant bits of the following number:
 *   x - 1 - 2^(W+B) + 2^W = x - 1 - 2^W (2^B - 1)
 *
 * Note that when W=0, the code is standard gamma-code.
 */
template <typename BitVector>
void WaveletTree<BitVector>::
fixedIntegerCode(BitVector& bits, uint32 x, uint32 w) {
  assert(x > 0);
  const uint64 wPow = static_cast<uint64>(1) << w;
  size_t B = floor(log(x - 1 + wPow)/log(2)) - w;
  //TODO: utils::logFloor
  for(size_t i = 0; i < B; ++i) bits.push_back(true);
  bits.push_back(false);

  uint64 y = x - fixedIntegerCodeTranslation(B,w);
  for(int i = w + B - 1; i >= 0; --i)
    bits.push_back( ((y >> i) & 1) == 1);
}

template <typename BitVector>
size_t WaveletTree<BitVector>::
fixedIntegerCodeTranslation(size_t bits, uint32 w) {
  return 1 + (((static_cast<size_t>(1) << bits) - 1) << w);
}

template <typename BitVector>
void WaveletTree<BitVector>::gammaCode(BitVector& bits, size_t integer) {
  int log = utils::logFloor(integer);
  for(int i = 0; i < log; ++i) {
    bits.push_back(true);
  }
  bits.push_back(false);
  for(int i = log-1; i >= 0; --i) {
    bits.push_back((integer >> i)&1);
  }
}

template <typename BitVector>
void WaveletTree<BitVector>::treeShape(BitVector& vec) const {
  size_t maxLen = 0;
  byte b = 0;
  std::vector<byte> symbols;
  for(size_t i = 0; i < 256; ++i, ++b) {
    size_t len = m_codes[b].size();
    if(len > 0) {
      symbols.push_back(b);
      if(len > maxLen) maxLen = len;
    }
  }
  // Largest symbol in alphabet
  utils::pushBits(vec, symbols.back(), 8);
  // Number of distinct symbols
  utils::pushBits(vec, symbols.size(), 8);

  int bytesInLongestCode;
  size_t packedInt = utils::packInteger(maxLen, &bytesInLongestCode);
  utils::pushBitsRev(vec, packedInt, 8*bytesInLongestCode);

  utils::binaryInterpolativeCode(symbols, symbols.back(), vec);

  for(size_t i = 0; i < symbols.size(); ++i) 
    utils::unaryCode(vec, maxLen - m_codes[symbols[i]].size() + 1);

#ifdef OPTIMIZED_INTEGER_CODE
  {
    maxLen = 0;
    std::vector<uint32> integers;
    for(typename std::map<uint32, BitVector>::const_iterator it = m_integerCodes.begin();
        it != m_integerCodes.end(); ++it)
    {
      integers.push_back(it->first);
      size_t len = it->second.size();
      if(len > maxLen) maxLen = len;
    }

    int bytes;
    // Largest symbol
    packedInt = utils::packInteger(integers.back(), &bytes);
    utils::pushBitsRev(vec, packedInt, 8*bytes);

    // Num of distinct symbols
    packedInt = utils::packInteger(integers.size(), &bytes);
    utils::pushBitsRev(vec, packedInt, 8*bytes);

    // Maximum length of code
    packedInt = utils::packInteger(maxLen, &bytes);
    utils::pushBitsRev(vec, packedInt, 8*bytes);
#ifndef SEMI_FIXED_CODE
    utils::binaryInterpolativeCode(integers, 0, integers.size() - 1, 1, integers.back(), vec);
#else    
    utils::binaryInterpolativeCode(integers, 0, integers.size() - 1, 0, integers.back(), vec);
#endif
    
    for(typename std::map<uint32, BitVector>::const_iterator it = m_integerCodes.begin();
        it != m_integerCodes.end(); ++it) {
      assert(it->second.size() > 0);
      utils::unaryCode(vec, maxLen - it->second.size() + 1);
    }

  }
#endif
}

template <typename BitVector>
#ifdef ENTROPY_PROFILER
template <typename Encoder, typename ProbabilisticModel, typename GammaModel,
          typename GapModel>
void WaveletTree<BitVector>::encodeTreeBF(Encoder& enc, ProbabilisticModel& pm,
                                          GammaModel& gm, GapModel& gapm)
#else
template <typename Encoder, typename ProbabilisticModel, typename GammaModel,
          typename GapModel>
void WaveletTree<BitVector>::encodeTreeBF(Encoder& enc, ProbabilisticModel& pm,
                                          GammaModel& gm, GapModel& gapm) const
#endif
{
  PROFILE("WaveletTree::encodeTreeBF");
  /* Additional bitvector is used to encode gaps and continuous runs in the
   * parent's bitvector. */
  typedef std::pair<TreeNode<BitVector>*, BitVector> InternalNode;
  /* TODO: If needing optimization helper bitvectors can be represented
   * using single larger bitvector. Now we just make redundant copies of
   * bitvectors.
   */

#ifdef ENTROPY_PROFILER
  uint64 bytesWritten = enc.counter();
#endif

  std::queue<InternalNode> queue;
  std::list<TreeNode<BitVector>*> integerCodeNodes;
  {
    // Root node
    InternalNode left, right;
    bool prev = !m_root->m_bitVector[0];
    for(size_t i = 0; i < m_root->m_bitVector.size(); ++i) {
      bool bit = m_root->m_bitVector[i];
      enc.encode(bit, pm.probabilityOfOne());
      pm.update(bit);
      BitVector& bv = bit? right.second : left.second;
      bv.push_back(prev != bit);
      prev = bit;
    }
    if(m_root->m_left) {
      if(m_root->m_left->m_hasSymbol) integerCodeNodes.push_back(m_root->m_left);
      else { left.first = m_root->m_left; queue.push(left); }
    }
    if(m_root->m_right) {
      if(m_root->m_right->m_hasSymbol) integerCodeNodes.push_back(m_root->m_right);
      else { right.first = m_root->m_right; queue.push(right); }
    }
  }


  // Nodes between root node and symbol nodes
  while(!queue.empty()) {
    pm.resetModel(); 
    gapm.resetModel();

    InternalNode left, right;
    InternalNode& node = queue.front();
    bool prev = !node.first->m_bitVector[0];
    if(node.first->m_left->m_hasSymbol || node.first->m_right->m_hasSymbol) {

      // Both children are symbol nodes, hence only "after gaps" are needed to encode
      if(node.first->m_left->m_hasSymbol && node.first->m_right->m_hasSymbol) {
        for(size_t i = 0; i  < node.first->m_bitVector.size(); ++i) {
          if(!node.second[i]) continue; //node.first->m_bitVector[i] is known
          bool bit = node.first->m_bitVector[i];

          enc.encode(bit, gapm.probabilityOfOne());
          gapm.update(bit);

        }
        integerCodeNodes.push_back(node.first->m_left);
        integerCodeNodes.push_back(node.first->m_right);
        
      } else if(node.first->m_left->m_hasSymbol) {
        right.first = node.first->m_right;

        for(size_t i = 0; i < node.first->m_bitVector.size(); ++i) {
          bool bit = node.first->m_bitVector[i];
          if(bit) right.second.push_back(prev != bit || node.second[i]);
          if(prev || node.second[i]) {
            if(node.second[i]) {
              enc.encode(bit, gapm.probabilityOfOne());
              gapm.update(bit);
              pm.updateState(bit);
            } else {
              enc.encode(bit, pm.probabilityOfOne());
              pm.update(bit);
            }
          }
          prev = bit; 
        }
        queue.push(right);
        integerCodeNodes.push_back(node.first->m_left);
      } else {
        // Shouldn't never happen, because we use canonical Huffman code
        assert(0);
      }
      
    } else {
      for(size_t i = 0; i < node.first->m_bitVector.size(); ++i) {
        bool bit = node.first->m_bitVector[i];
        if(node.second[i]) {
          enc.encode(bit, gapm.probabilityOfOne());
          gapm.update(bit);
          pm.updateState(bit);
        } else {
          enc.encode(bit, pm.probabilityOfOne());
          pm.update(bit);
        }
        BitVector& bv = bit? right.second : left.second;
        bv.push_back(prev != bit || node.second[i]);
        prev = bit;
      }
      left.first = node.first->m_left;
      queue.push(left);
      right.first = node.first->m_right;
      queue.push(right);
    }
    queue.pop();
  }

#ifdef ENTROPY_PROFILER
  uint64 nBytes = enc.counter();
  m_bytesForCharacters += (nBytes - bytesWritten);
  bytesWritten = nBytes;
#endif
  
  // Synchronized integer-coding phase (symbol nodes and integer-code nodes)
  std::list<TreeNode<BitVector>*> left, right;
  while(!integerCodeNodes.empty() || !left.empty() || !right.empty()) {
    integerCodeNodes.splice(integerCodeNodes.end(), left);
    integerCodeNodes.splice(integerCodeNodes.end(), right);
    gm.resetModel();
    while(!integerCodeNodes.empty()) {

      TreeNode<BitVector>* node = integerCodeNodes.front();
      for(size_t i = 0; i < node->m_bitVector.size(); ++i) {
        enc.encode(node->m_bitVector[i], gm.probabilityOfOne());
        gm.update(node->m_bitVector[i]);
      }
      integerCodeNodes.pop_front();
#ifdef OPTIMIZED_INTEGER_CODE
#ifndef SEMI_FIXED_CODE      
      if(node->m_left && !node->m_left->m_hasSymbol)
        left.push_back(node->m_left);
      if(node->m_right && !node->m_right->m_hasSymbol)
        right.push_back(node->m_right);
#else
      if(node->m_left &&
         (!node->m_left->m_hasSymbol ||
          (node->m_left->m_hasSymbol && node->m_left->m_symbol == 0)))
        left.push_back(node->m_left);
      if(node->m_right &&
         (!node->m_right->m_hasSymbol ||
          (node->m_right->m_hasSymbol && node->m_right->m_symbol == 0)))
        right.push_back(node->m_right);
#endif      
#else      
      if(node->m_left)
        left.push_back(node->m_left);
      if(node->m_right)
        right.push_back(node->m_right);
#endif
    }
  }
#ifdef ENTROPY_PROFILER
  nBytes = enc.counter();
  m_bytesForRuns += (nBytes - bytesWritten);
#endif
}

template <typename BitVector>
template <typename Encoder, typename ProbabilisticModel>
void WaveletTree<BitVector>::encodeTree(Encoder& enc, ProbabilisticModel& pm) const
{
  encodeTree(m_root, enc, pm);
}

template <typename BitVector>
template <typename Encoder, typename ProbabilisticModel>
void WaveletTree<BitVector>::encodeTree(TreeNode<BitVector> *node,
                                        Encoder& enc,
                                        ProbabilisticModel& pm) const
{
  for(size_t i = 0; i < node->m_bitVector.size(); ++i) {
    enc.encode(node->m_bitVector[i], pm.probabilityOfOne());
    pm.update(node->m_bitVector[i]);
  }
  if(node->m_left) encodeTree(node->m_left, enc, pm);
  if(node->m_right) encodeTree(node->m_right, enc, pm);
}

/** We have to detect the end of each gamma-code by decoding that code.
 *  So to decode the whole tree two additional decoding-methods for gamma
 *  codes are needed. This one is used when the first 0 is encountered
 *  (we have to keep track how many bits in code to read). 
 */
template <typename BitVector>
template <typename Decoder, typename ProbabilisticModel>
void WaveletTree<BitVector>::decodeTreeGammaDec(TreeNode<BitVector>* node,
                                                size_t nodeSize,
                                                size_t bitsInCode,
                                                Decoder& dec,
                                                ProbabilisticModel& pm) const
{
  if(bitsInCode == 0) return;
  size_t ones = 0;
  for(size_t i = 0; i < nodeSize; ++i) {
    bool bit = dec.decode(pm.probabilityOfOne());
    pm.update(bit);
    if(bit) ++ones;
    node->m_bitVector.push_back(bit);
  }
  if(nodeSize - ones > 0) {
    if(!node->m_left) node->m_left = new TreeNode<BitVector>();
    decodeTreeGammaDec(node->m_left, nodeSize - ones, bitsInCode-1, dec, pm);
  }
  if(ones > 0) {
    if(!node->m_right) node->m_right = new TreeNode<BitVector>();
    decodeTreeGammaDec(node->m_right, ones, bitsInCode-1, dec, pm);
  }
}

/** This is used when we haven't yet seen the first 0-bit in gamma-code.
 */
template <typename BitVector>
template <typename Decoder, typename ProbabilisticModel>
void WaveletTree<BitVector>::decodeTreeGammaInc(TreeNode<BitVector>* node,
                                                size_t nodeSize,
                                                size_t bitsInCode,
                                                Decoder& dec,
                                                ProbabilisticModel& pm) const
{
  size_t ones = 0;
  for(size_t i = 0; i < nodeSize; ++i) {
    bool bit = dec.decode(pm.probabilityOfOne());
    pm.update(bit);
    if(bit) ++ones;
    node->m_bitVector.push_back(bit);
  }
  if(nodeSize - ones > 0) {
    if(!node->m_left) node->m_left = new TreeNode<BitVector>();
    decodeTreeGammaDec(node->m_left, nodeSize - ones, bitsInCode, dec, pm);
  }
  if(ones > 0) {
    if(!node->m_right) node->m_right = new TreeNode<BitVector>();
    decodeTreeGammaInc(node->m_right, ones, bitsInCode +1, dec, pm);
  }
}

/**Used in decoding the tree. */
template<typename BitVector>
#ifndef OPTIMIZED_INTEGER_CODE
// Used for gamma codes, tracks progression of the code with variables
struct IntegerNode {
  IntegerNode(TreeNode<BitVector>* node, size_t bits, size_t len, byte status)
      : m_node(node), m_bits(bits), m_gammaLen(len), m_gammaStatus(status) {}
  TreeNode<BitVector> *m_node;
  size_t m_bits;
  size_t m_gammaLen;
  byte m_gammaStatus; // Flag telling in what phase of gamma-code decoding is
};
#else
#ifndef SEMI_FIXED_CODE
// Used for more general integer codes. Tracks progression of the code by
// following the corresponding integer code tree
struct IntegerNode {
  IntegerNode(TreeNode<BitVector>* node, TreeNode<BitVector>* intNode, size_t bits)
      : m_node(node), m_intNode(intNode), m_bits(bits) {}
  TreeNode<BitVector> *m_node;
  TreeNode<BitVector> *m_intNode;
  size_t m_bits;
};
#else
// Used for semi-fixed codes. Here we need to be able to point where we are in
// integer code tree. In addition it is needed to know if we are in the fixed-code
// part and how many leading ones the code has.
struct IntegerNode {
  IntegerNode(TreeNode<BitVector>* node, TreeNode<BitVector>* intNode, size_t bits,
              uint32 leadingOnes, byte codeStatus)
      : m_node(node), m_intNode(intNode), m_bits(bits), m_leadingOnes(leadingOnes),
        m_codeStatus(codeStatus) {}
  TreeNode<BitVector> *m_node;
  TreeNode<BitVector> *m_intNode;
  size_t m_bits;
  uint32 m_leadingOnes;
  byte m_codeStatus;
};
#endif
#endif

/** Note that we assume that the wavelet tree has canonical Huffman shape ie.
 *  the longest codewords are placed to right.
 */
template<typename BitVector>
template <typename Decoder, typename ProbabilisticModel, typename GammaModel,
          typename GapModel>
void WaveletTree<BitVector>::decodeTreeBF(size_t rootSize,
                                          Decoder& dec,
                                          ProbabilisticModel& pm,
                                          GammaModel& gm,
                                          GapModel& gapm) const
{
  //TODO: If needed optimize redundant copying of bitvectors!
  typedef std::pair<TreeNode<BitVector>*, BitVector> InternalNode; 
  
  std::queue<InternalNode> queue;
  std::list<IntegerNode<BitVector> > integerCodeNodes;

  //Decoding of the root node
  {
    BitVector left, right;

    bool prev = dec.decode(pm.probabilityOfOne());
    pm.update(prev);
    m_root->m_bitVector.push_back(prev);

    if(prev) right.push_back(true);
    else left.push_back(true);

    for(size_t i = 1; i < rootSize; ++i) {
      bool bit = dec.decode(pm.probabilityOfOne());
      pm.update(bit);
      m_root->m_bitVector.push_back(bit);
      BitVector& gVector = bit? right: left;
      gVector.push_back(prev != bit);
      prev = bit;
    }
    // Left node has to always exist
    if(m_root->m_left->m_hasSymbol) {
      integerCodeNodes.push_back(
#ifndef OPTIMIZED_INTEGER_CODE
        IntegerNode<BitVector>(m_root->m_left, left.size(), 0, 0));
#else
#ifndef SEMI_FIXED_CODE
        IntegerNode<BitVector>(m_root->m_left, m_integerCodeTree, left.size()));
#else
        IntegerNode<BitVector>(m_root->m_left, m_integerCodeTree, left.size(),
                               0, 0));
#endif
#endif      
    } else {
      queue.push(std::make_pair(m_root->m_left, left));
    }
    
    if(right.size() > 0) {
      if(m_root->m_right->m_hasSymbol) {
        integerCodeNodes.push_back(
#ifndef OPTIMIZED_INTEGER_CODE
          IntegerNode<BitVector>(m_root->m_right, right.size(), 0, 0));
#else
#ifndef SEMI_FIXED_CODE
          IntegerNode<BitVector>(m_root->m_right, m_integerCodeTree, right.size()));
#else
          IntegerNode<BitVector>(m_root->m_right, m_integerCodeTree, right.size(),
                                 0, 0));
#endif
#endif      
      } else {
        queue.push(std::make_pair(m_root->m_right, right));
      }
    }
  }

  //Decoding of the internal nodes
  {
    while(!queue.empty()) {
      pm.resetModel();
      gapm.resetModel();

      BitVector left, right;
      InternalNode& node = queue.front();
      // Node must have both left and right child
      if(node.first->m_left->m_hasSymbol || node.first->m_right->m_hasSymbol) {
        if(node.first->m_left->m_hasSymbol && node.first->m_right->m_hasSymbol) {
          size_t ones = 0;
          bool prev = true;
          for(size_t i = 0; i < node.second.size(); ++i) {
            if(!node.second[i]) {
              prev = !prev;
            } else {
              prev = dec.decode(gapm.probabilityOfOne());
              gapm.update(prev);
            }
            node.first->m_bitVector.push_back(prev);
            if(prev) ++ones;
          }
          integerCodeNodes.push_back(IntegerNode<BitVector>(
#ifndef OPTIMIZED_INTEGER_CODE
              node.first->m_left, node.second.size() - ones, 0, 0));
#else
#ifndef SEMI_FIXED_CODE
              node.first->m_left, m_integerCodeTree, node.second.size() - ones));
#else
      node.first->m_left, m_integerCodeTree, node.second.size() - ones,0,0));
#endif
#endif
          integerCodeNodes.push_back(IntegerNode<BitVector>(
#ifndef OPTIMIZED_INTEGER_CODE
              node.first->m_right, ones, 0, 0));
#else
#ifndef SEMI_FIXED_CODE
              node.first->m_right, m_integerCodeTree, ones));
#else
              node.first->m_right, m_integerCodeTree, ones, 0, 0));
#endif
#endif
        } else {
          assert(!node.first->m_right->m_hasSymbol);
          bool prev = true;
          for(size_t i = 0; i < node.second.size(); ++i) {
            bool bit;
            if(!node.second[i] && !prev) {
              bit = true;
            } else if(node.second[i]) {
              bit = dec.decode(gapm.probabilityOfOne());
              gapm.update(bit);
              pm.updateState(bit);
            } else {
              bit = dec.decode(pm.probabilityOfOne());
              pm.update(bit);
            }
            node.first->m_bitVector.push_back(bit);
            if(bit) right.push_back(prev != bit || node.second[i]);
            prev = bit;
          }
          integerCodeNodes.push_back(IntegerNode<BitVector>(
#ifndef OPTIMIZED_INTEGER_CODE
            node.first->m_left, node.second.size() - right.size(), 0, 0));
#else
#ifndef SEMI_FIXED_CODE
            node.first->m_left, m_integerCodeTree, node.second.size() - right.size()));
#else
            node.first->m_left, m_integerCodeTree, node.second.size() - right.size(),0 ,0));
#endif
#endif
          queue.push(std::make_pair(node.first->m_right, right));
        }
      } else { //both children are also internal nodes
        bool prev = true;
        for(size_t i = 0; i < node.second.size(); ++i) {
          bool bit;
          if(node.second[i]) {
            bit = dec.decode(gapm.probabilityOfOne());
            gapm.update(bit);
            pm.updateState(bit);
          } else {
            bit = dec.decode(pm.probabilityOfOne());
            pm.update(bit);
          }
          node.first->m_bitVector.push_back(bit);
          BitVector& gapVector = bit? right: left;
          gapVector.push_back(prev != bit || node.second[i]);
          prev = bit;
        }
        queue.push(std::make_pair(node.first->m_left, left));
        queue.push(std::make_pair(node.first->m_right, right));
      }
      queue.pop();
    }
  }

  //Decoding of integer-code nodes
  {
    std::list<IntegerNode<BitVector> > left, right;
    while(!integerCodeNodes.empty() || !left.empty() || !right.empty()) {
      gm.resetModel();
      
      while(!integerCodeNodes.empty()) {
        IntegerNode<BitVector> node = integerCodeNodes.front();
        integerCodeNodes.pop_front();
#ifdef OPTIMIZED_INTEGER_CODE
#ifndef SEMI_FIXED_CODE
        assert(node.m_intNode);
        if(node.m_intNode->m_hasSymbol) {
          node.m_node->m_hasSymbol = true;
          node.m_node->m_symbol = node.m_intNode->m_symbol;
          continue;
        }
#else
        if(node.m_intNode && node.m_intNode->m_hasSymbol) {
          node.m_node->m_hasSymbol = true;
          node.m_node->m_symbol = node.m_intNode->m_symbol;
          if(node.m_node->m_symbol != 0) continue;
        }
#endif
#endif        
        size_t ones = 0;
        for(size_t i = 0; i < node.m_bits; ++i) {
          bool bit = dec.decode(gm.probabilityOfOne());
          gm.update(bit);
          if(bit) ++ones;
          node.m_node->m_bitVector.push_back(bit);
        }
        
#ifndef OPTIMIZED_INTEGER_CODE
        if(node.m_gammaStatus == 2 && node.m_gammaLen == 1) continue;

        
        if(node.m_bits > ones && node.m_gammaStatus != 0) {
          if(!node.m_node->m_left)
            node.m_node->m_left = new TreeNode<BitVector>();

          IntegerNode<BitVector> lnode(node.m_node->m_left, node.m_bits - ones,
                                     node.m_gammaLen, node.m_gammaStatus);
          if(node.m_gammaStatus == 1) {
            lnode.m_gammaStatus = 2;
          } else if(node.m_gammaStatus == 2) {
            --lnode.m_gammaLen;
          }
          left.push_back(lnode);
        }

        if(ones > 0) {
          if(!node.m_node->m_right)
            node.m_node->m_right = new TreeNode<BitVector>();

          IntegerNode<BitVector> rnode(node.m_node->m_right, ones,
                                       node.m_gammaLen, node.m_gammaStatus);
          if(node.m_gammaStatus == 0) {
            rnode.m_gammaLen = 1;
            rnode.m_gammaStatus = 1;
          } else if(node.m_gammaStatus == 1) {
            ++rnode.m_gammaLen;
          } else if(node.m_gammaStatus == 2) {
            --rnode.m_gammaLen;
          }
          right.push_back(rnode);
        }
#else
        if(node.m_bits > ones) {
          if(!node.m_node->m_left)
            node.m_node->m_left = new TreeNode<BitVector>();

#ifndef SEMI_FIXED_CODE
          assert(node.m_intNode->m_left);
          IntegerNode<BitVector> lnode(node.m_node->m_left,
                                       node.m_intNode->m_left,
                                       node.m_bits - ones);
          left.push_back(lnode);
#else

          IntegerNode<BitVector> lnode(node.m_node->m_left,
                                       node.m_intNode?node.m_intNode->m_left:0,
                                       node.m_bits - ones, node.m_leadingOnes,
                                       node.m_codeStatus);

          if(!node.m_intNode || !node.m_intNode->m_left) {
            if(node.m_codeStatus == 0) {
              lnode.m_codeStatus = 2;
              lnode.m_leadingOnes = m_W;
            } else if(node.m_codeStatus == 1) {
              lnode.m_codeStatus = 2;
              lnode.m_leadingOnes += m_W;
            } else {
              --lnode.m_leadingOnes;
            }
          }

          if(lnode.m_codeStatus != 2 || lnode.m_leadingOnes > 0) {
            left.push_back(lnode);
          }
          
#endif
        }

        if(ones > 0) {
          if(!node.m_node->m_right)
            node.m_node->m_right = new TreeNode<BitVector>();

#ifndef SEMI_FIXED_CODE
          assert(node.m_intNode->m_right);
          IntegerNode<BitVector> rnode(node.m_node->m_right,
                                       node.m_intNode->m_right,
                                       ones);
          right.push_back(rnode);
#else

          IntegerNode<BitVector> rnode(node.m_node->m_right,
                                       node.m_intNode?node.m_intNode->m_right:0,
                                       ones, node.m_leadingOnes,
                                       node.m_codeStatus);

          if(!node.m_intNode || !node.m_intNode->m_right) {
            if(node.m_codeStatus == 0) {
              rnode.m_codeStatus = 1;
              ++rnode.m_leadingOnes;
            } else if(node.m_codeStatus == 1) {
              ++rnode.m_leadingOnes;
            } else {
              --rnode.m_leadingOnes;
            }
          }
          
          if(rnode.m_codeStatus != 2 || rnode.m_leadingOnes > 0) {
            right.push_back(rnode);
          }
#endif
        }
#endif        
      }
      integerCodeNodes.splice(integerCodeNodes.end(), left);
      integerCodeNodes.splice(integerCodeNodes.end(), right);
    }
  }
}

template<typename BitVector>
template <typename Decoder, typename ProbabilisticModel>
void WaveletTree<BitVector>::decodeTree(size_t rootSize,
                                        Decoder& dec,
                                        ProbabilisticModel& pm) const
{
  decodeTree(m_root, rootSize, dec, pm);
}

template <typename BitVector>
template <typename Decoder, typename ProbabilisticModel>
void WaveletTree<BitVector>::decodeTree(TreeNode<BitVector>* node, size_t nodeSize,
                                        Decoder& dec, ProbabilisticModel& pm) const
{
  if(!node->m_hasSymbol) {
    size_t ones = 0;
    for(size_t i = 0; i < nodeSize; ++i) {
      bool bit = dec.decode(pm.probabilityOfOne());
      pm.update(bit);
      if(bit) ++ones;
      node->m_bitVector.push_back(bit);
    }
    if(nodeSize - ones > 0) decodeTree(node->m_left, nodeSize - ones, dec, pm);
    if(ones > 0) decodeTree(node->m_right, ones, dec, pm);
  } else {
    decodeTreeGammaInc(node, nodeSize, 0, dec, pm);
  }
}
  
/** Pushes bitvector to tree by starting from the root. Assumes that every
 *  node on the path exists. This is used when pushing characters into the tree
 *  whose shape is already formed and codes are chosen by this shape.
 *
 * @param bits BitVector representing the symbol to be pushed.
 * @return Node which is at the depth bits.size() when the root is at depth 0.
 *         This node should have symbol (ie. m_hasSymbol == true).                   
 */
template <typename BitVector>
TreeNode<BitVector> *WaveletTree<BitVector>::pushBits(const BitVector& bits)
{
  TreeNode<BitVector> *node = m_root;
  for(size_t i = 0; i < bits.size(); ++i) {
    node->m_bitVector.push_back(bits[i]);
    if(bits[i]) {
      assert(node->m_right);
      node = node->m_right;
    } else {
      assert(node->m_left);
      node = node->m_left;
    }
  }
  return node;
}


/** Pushes bitvector to tree by starting the node given. Creates new nodes if
 *  they don't exist on the path. However this doesn't create the leaf node which
 *  would come after storing all of the bits in their respectable nodes.
 *
 * @param node Node to start.
 * @param bits Bits to push to the tree.
 */
template <typename BitVector>
void
WaveletTree<BitVector>::pushBits(TreeNode<BitVector> *node, const BitVector& bits)
{
  assert(node);
  assert(bits.size() > 0);
  for(size_t i = 0; i < bits.size() - 1; ++i) {
    node->m_bitVector.push_back(bits[i]);
    if(bits[i]) {
      if(!node->m_right) node->m_right = new TreeNode<BitVector>();
      node = node->m_right;
    } else {
      if(!node->m_left) node->m_left = new TreeNode<BitVector>();
      node = node->m_left;
    }
  }
  node->m_bitVector.push_back(bits.back());
}

template <typename BitVector> void
WaveletTree<BitVector>::pushBits(TreeNode<BitVector> *node, const BitVector& bits,
                                 uint32 symbol)
{
  assert(node);
  assert(bits.size() > 0);
  for(size_t i = 0; i < bits.size() - 1; ++i) {
    node->m_bitVector.push_back(bits[i]);
    if(bits[i]) {
      if(!node->m_right) node->m_right = new TreeNode<BitVector>();
      node = node->m_right;
    } else {
      if(!node->m_left) node->m_left = new TreeNode<BitVector>();
      node = node->m_left;
    }
  }
  node->m_bitVector.push_back(bits.back());
  if(bits.back() && !node->m_right)
    node->m_right = new TreeNode<BitVector>(symbol);
  else if(!bits.back() && !node->m_left)
    node->m_left = new TreeNode<BitVector>(symbol);
}

template <typename BitVector>
void WaveletTree<BitVector>::pushRun(byte symbol, size_t runLength)
{
  TreeNode<BitVector> *node = pushBits(m_codes[symbol]);
  assert(node->m_hasSymbol);
#ifndef OPTIMIZED_INTEGER_CODE
  BitVector integerCode;
  gammaCode(integerCode, runLength);
  pushBits(node, integerCode);
#else

#ifdef SEMI_FIXED_CODE
  if(m_integerCodes.find(runLength) == m_integerCodes.end()) {
    BitVector integerCode(m_integerCodes[0]);
    fixedIntegerCode(integerCode, runLength, m_W);
    pushBits(node, integerCode, runLength);
  } else {
    pushBits(node, m_integerCodes[runLength], runLength);
  }
#else
  pushBits(node, m_integerCodes[runLength], runLength);
#endif

#endif  
}

template <typename BitVector> template <typename OutputIterator>
size_t WaveletTree<BitVector>::message(OutputIterator out) const {
  size_t len = 0;
  size_t msgSize = m_root->m_bitVector.size();
  std::map<TreeNode<BitVector>*, size_t> bitsSeen;
  for(size_t j = 0; j < msgSize; ++j) {
    bool bit;
    TreeNode<BitVector> *node = m_root;
    size_t i = j;
    do {
      bit = node->m_bitVector[i];
      // Update bits seen and perform rank
      if(bit) i = bitsSeen[node]++;
      else i = i - bitsSeen[node]; 
      node = bit?node->m_right:node->m_left;
    } while(!node->m_hasSymbol);
    byte symbol = node->m_symbol;

    // Decoding of integer-code
#ifndef OPTIMIZED_INTEGER_CODE
    size_t runBits = 0;
    bit = node->m_bitVector[i];
    while(bit) {
      i = bitsSeen[node]++;
      node = bit?node->m_right:node->m_left;
      ++runBits;
      bit = node->m_bitVector[i];
    }
    size_t runLength = 1;
    for(size_t k = 0; k < runBits; ++k) {
      runLength <<= 1;
      if(bit) i = bitsSeen[node]++;
      else i = i - bitsSeen[node]; 
      node = bit?node->m_right:node->m_left;
      bit = node->m_bitVector[i];
      runLength |= (bit?1:0);
    }
#else
    do {
      bit = node->m_bitVector[i];

      if(bit) {
        i = bitsSeen[node]++;
        assert(node->m_right);
        node = node->m_right;
      } else {
        i = i - bitsSeen[node];
        assert(node->m_left);
        node = node->m_left;
      }
    } while(!node->m_hasSymbol);
    size_t runLength = node->m_symbol;

#ifdef SEMI_FIXED_CODE
    size_t leadingOnes = 0;
    if(runLength == 0) {
      bit = node->m_bitVector[i];
      while(bit) {
        ++leadingOnes;
        i = bitsSeen[node]++;
        node = bit?node->m_right:node->m_left;
        bit = node->m_bitVector[i];
      }
      
      for(size_t k = 0; k < leadingOnes + m_W; ++k) {
        runLength <<= 1;
        if(bit) i = bitsSeen[node]++;
        else i = i - bitsSeen[node];
        node = bit?node->m_right:node->m_left;
        bit = node->m_bitVector[i];
        runLength |= (bit?1:0);
      }
      runLength += fixedIntegerCodeTranslation(leadingOnes, m_W);
    }
#endif

#endif //OPTIMIZED_INTEGER_CODE
    for(size_t k = 0; k < runLength; ++k) {
      *out++ = symbol;
    }
    len += runLength;
  }
  return len;
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
void WaveletTree<BitVector>::assignPrefixCodes(
    std::vector<std::pair<uint64, uint32> >& lengths)
{
  std::sort(lengths.begin(), lengths.end());
  m_root = new TreeNode<BitVector>();
  assignPrefixCodes(lengths, m_root, 0, 0);
}

template <typename BitVector>
size_t WaveletTree<BitVector>::assignPrefixCodes(
    std::vector<std::pair<uint64, uint32> >& lengths, TreeNode<BitVector>* node,
    size_t elem, size_t bits)
{
  if(elem >= lengths.size()) return elem;
  /*if(bits == lengths[elem].first) {
    node->m_hasSymbol = true;
    node->m_symbol = lengths[elem].second;
    return elem + 1;
  } else */
  if(bits == lengths[elem].first - 1) {
    if(!node->m_left) {
      node->m_left = new TreeNode<BitVector>(lengths[elem].second);
      elem = assignPrefixCodes(lengths, node, elem+1, bits);
    } else {
      assert(!node->m_right);
      node->m_right = new TreeNode<BitVector>(lengths[elem].second);
      ++elem;
    }
    return elem;
  }
  assert(bits < lengths[elem].first - 1);
  if(!node->m_left) {
    node->m_left = new TreeNode<BitVector>();
    elem = assignPrefixCodes(lengths, node->m_left, elem, bits + 1);
  }
  assert(!node->m_right);
  if(elem < lengths.size()) {
    node->m_right = new TreeNode<BitVector>();
    elem = assignPrefixCodes(lengths, node->m_right, elem, bits + 1);
  }
  return elem;
}

template <typename BitVector>
void WaveletTree<BitVector>::pushMessage(const byte* src, size_t length) {
  const byte *prev = src;
  const byte *curr = src+1;
  do{
    while(curr < src + length && *prev == *curr) ++curr;
    pushRun(*prev, curr-prev);
    prev = curr++;
  } while(curr < src + length);
  if(prev < src + length) pushRun(*prev,1);
}

template <typename BitVector>
TreeNode<BitVector>*
WaveletTree<BitVector>::createHuffmanShape(const uint64 *runFreqs)
{
  MinimumHeap<TreeNode<BitVector>* > heap(256);
  byte b = 0U;
  for(int i = 0; i < 256; ++i,++b) {
    if(runFreqs[b] != 0) {
      TreeNode<BitVector> *node = new TreeNode<BitVector>(b);
      heap.insert(node, runFreqs[b]);
    }
  }
  if(heap.size() == 1) {
    TreeNode<BitVector> *root = new TreeNode<BitVector>(heap.deleteMin().first, 0);
    heap.insert(root, 1);
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
template <typename BVectors>
void WaveletTree<BitVector>::collectCodes(BVectors& codes, TreeNode<BitVector> *root)
{
  BitVector vec;
  collectCodes(codes, vec, root);
}

template <typename BitVector>
template <typename BVectors> void
WaveletTree<BitVector>::collectCodes(BVectors& codes, BitVector& vec,
                                     TreeNode<BitVector> *node)
{
  if(node->m_left == 0 && node->m_right == 0) {
    codes[node->m_symbol] = vec;
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

template <typename BitVector>
uint32 WaveletTree<BitVector>::findParametersForSemiFixedCodes(
    std::vector<std::pair<uint64, uint32> >& integerFrequencies,
    size_t totalFrequencies) {

    std::sort(integerFrequencies.rbegin(), integerFrequencies.rend());

    size_t rejectedSum = 0;
    const size_t common = 10;
    while(true) {
      if(integerFrequencies.back().first > common || integerFrequencies.size() <= 1) break;
      rejectedSum += integerFrequencies.back().first;
      integerFrequencies.pop_back();
    }
    std::reverse(integerFrequencies.begin(), integerFrequencies.end());
    // '0' Is used for the symbols which have semifixed code
    if(rejectedSum > 0) {
      integerFrequencies.push_back(std::make_pair(rejectedSum, 0));
      int i = integerFrequencies.size() - 2;
      while(i >= 0 && rejectedSum < integerFrequencies[i].first) {
        std::swap(integerFrequencies[i], integerFrequencies[i+1]);
        --i;
      }
    }

    return m_W;
}


} //namespace bwtc
#endif
