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
#include <cassert>
#include <map>
#include <queue>
#include <utility>
#include <vector>

namespace bwtc {

template <typename BitVector>
struct TreeNode {
  TreeNode() : m_left(0), m_right(0), m_hasSymbol(false) {}
  TreeNode(byte symbol) : m_left(0), m_right(0), m_hasSymbol(true),
                          m_symbol(symbol) {}
  TreeNode(TreeNode<BitVector>* left, TreeNode<BitVector>* right)
      : m_left(left), m_right(right), m_hasSymbol(false) {}
  ~TreeNode() {}

  size_t rank(bool bit, size_t i) const;
  size_t totalBits() const;

  BitVector m_bitVector;
  TreeNode<BitVector>* m_left;
  TreeNode<BitVector>* m_right;
  bool m_hasSymbol;
  byte m_symbol;
  
};

/**If ever using something suitable for fast rank-queries, there should be
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

  /**Writes the encoding of the shape of tree into given vector.*/
  void treeShape(BitVector& vector) const;

  /**Writes encoding of the shape of tree into given vector. The shape of
   * the tree is in form presented in "Housekeeping for Prefix Codes" by
   * Turpin and Moffat. */
  void treeShape2(BitVector& vector) const;

  size_t bitsInRoot() const { return m_root->m_bitVector.size(); }
  size_t totalBits() const { return m_root->totalBits(); }

  template <typename OutputIterator>
  void message(OutputIterator out) const;

  /**Input is required to have readBit()-method returning bool.*/
  template<typename Input>
  size_t readShape(Input& input);

  /**Input is required to have readBit()-method returning bool and readByte().
   */
  template<typename Input>
  size_t readShape2(Input& input);

  /** Encoder is required to have encode(bool bit, Probability prob)-method.
    * ProbabilisticModel is required to have probabilityOfOne()- and
    * update(bool bit)-methods. */
  template <typename Encoder, typename ProbabilisticModel>
  void encodeTree(Encoder& enc, ProbabilisticModel& pm) const;

  template <typename Encoder, typename ProbabilisticModel>
  void encodeTreeBF(Encoder& enc, ProbabilisticModel& pm) const;

  /** Decoder is required to have decode(Probability prob)-method.
    * ProbabilisticModel is required to have probabilityOfOne()- and
    * update(bool bit)-methods. */
  template <typename Decoder, typename ProbabilisticModel>
  void decodeTree(size_t rootSize, Decoder& enc, ProbabilisticModel& pm) const;

  template <typename Decoder, typename ProbabilisticModel>
  void decodeTreeBF(size_t rootSize, Decoder& enc, ProbabilisticModel& pm) const;

  const BitVector& code(byte symbol) { return m_codes[symbol]; }
  
  static TreeNode<BitVector> *createHuffmanShape(const uint64 *runFreqs);
  static void collectCodes(BitVector *codes, TreeNode<BitVector> *root);
  static void collectCodes(BitVector *codes, BitVector& vec,
                           TreeNode<BitVector> *node);
  static void pushBits(TreeNode<BitVector> *node, const BitVector& bits);
  static void gammaCode(BitVector& bits, size_t integer);

  /**Assigns prefix codes based on the given lengths. Actual codes are stored
   * into m_codes.
   *
   * @param lengths Lengths of the codewords to be assigned given as
   *                <length, symbol>-pairs. This array WILL be modified.
   */
  void assignPrefixCodes(std::vector<std::pair<uint64, byte> >& lengths);

  
 private:
  TreeNode<BitVector>* m_root;
  BitVector m_codes[256];

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

  template <typename Input>
  TreeNode<BitVector> *readShapeDfs(Input& input, const std::vector<byte> symbols,
                                    size_t *bitsUsed);

  void outputShapeDfs(BitVector& output, size_t depth,
                      const std::vector<byte>& symbols) const;

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
  size_t assignPrefixCodes(std::vector<std::pair<uint64, byte> >& lengths,
                         TreeNode<BitVector>* node, size_t elem, size_t bits);


};

template <typename BitVector>
WaveletTree<BitVector>::WaveletTree() {
  m_root = new TreeNode<BitVector>();
}

template <typename BitVector>
WaveletTree<BitVector>::WaveletTree(const byte *src, size_t length) {
  uint64 runFreqs[256] = {0};
  utils::calculateRunFrequencies(runFreqs, src, length);

#if 1
  std::vector<std::pair<uint64, byte> > codeLengths;
  utils::calculateHuffmanLengths(codeLengths, runFreqs);

  assignPrefixCodes(codeLengths);
  collectCodes(m_codes, m_root);

#else
  // Old way to construct huffman shape
  m_root = createHuffmanShape(runFreqs);
  collectCodes(m_codes, m_root);
#endif

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
WaveletTree<BitVector>::~WaveletTree() {
  destroy(m_root);
}

template <typename BitVector>
void WaveletTree<BitVector>::destroy(TreeNode<BitVector>* node) {
  if(node->m_left) destroy(node->m_left);
  if(node->m_right) destroy(node->m_right);
  delete node;
}


template <typename BitVector> template <typename Input>
size_t WaveletTree<BitVector>::readShape2(Input& input) {
  assert(m_root->m_left == 0 && m_root->m_right == 0);
  size_t maxSym = input.readByte();
  size_t symbols = input.readByte();
  size_t bitsRead = 16;
  size_t maxLen = 0;
  size_t read = 0xff;
  size_t j = 0;
  while(read & 0x80) {
    read = input.readByte();
    maxLen |= ((read & 0x7f) << j);
    j += 7;
    bitsRead += 8;
  }
  std::vector<byte> alphabet;
  bitsRead += utils::binaryInterpolativeDecode(alphabet, input,
                                               maxSym, symbols);
  std::vector<std::pair<uint64, byte> > codeLengths;
  for(size_t i = 0; i < symbols; ++i) {
    size_t n = utils::unaryDecode(input);
    bitsRead += n;
    size_t len = maxLen - n + 1;
    codeLengths.push_back(std::make_pair(len, alphabet[i]));
  }

  assignPrefixCodes(codeLengths);
  collectCodes(m_codes, m_root);
  
  return bitsRead;
}

template <typename BitVector> template <typename Input>
size_t WaveletTree<BitVector>::readShape(Input& input) {
  assert(m_root->m_left == 0 && m_root->m_right == 0);
  std::vector<byte> symbols;
  byte b = 0;
  for(size_t i = 0; i < 256; ++i, ++b) {
    if(input.readBit()) symbols.push_back(b);
  }
  std::vector<byte> leftSymbols, rightSymbols;
  for(size_t i = 0; i < symbols.size(); ++i) {
    bool bit = input.readBit();
    m_codes[symbols[i]].push_back(bit);

    if(bit) rightSymbols.push_back(symbols[i]);
    else leftSymbols.push_back(symbols[i]);
  }
  size_t bitsRead = 256 + symbols.size();
  if(leftSymbols.size() > 0) {
    m_root->m_left = readShapeDfs(input, leftSymbols, &bitsRead);
  }
  if(rightSymbols.size() > 0) {
    m_root->m_right = readShapeDfs(input, rightSymbols, &bitsRead);
  }
  return bitsRead;
}

template <typename BitVector> template <typename Input>
TreeNode<BitVector>*
WaveletTree<BitVector>::readShapeDfs(Input& input, const std::vector<byte> symbols,
                                     size_t *bitsUsed)
{
  if(symbols.size() == 1) {
    return new TreeNode<BitVector>(symbols[0]);
  }
  *bitsUsed += symbols.size();
  TreeNode<BitVector> *node = new TreeNode<BitVector>();
  std::vector<byte> leftSymbols, rightSymbols;
  for(size_t i = 0; i < symbols.size(); ++i) {
    bool bit = input.readBit();
    m_codes[symbols[i]].push_back(bit);

    if(bit) rightSymbols.push_back(symbols[i]);
    else leftSymbols.push_back(symbols[i]);
  }
  if(leftSymbols.size() > 0) {
    node->m_left = readShapeDfs(input, leftSymbols, bitsUsed);
  }
  if(rightSymbols.size() > 0) {
    node->m_right = readShapeDfs(input, rightSymbols, bitsUsed);
  }
  return node;
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
void WaveletTree<BitVector>::treeShape(BitVector& vector) const {
  byte b = 0;
  std::vector<byte> symbols;
  for(size_t i = 0; i < 256; ++i, ++b) {
    if(m_codes[b].size() > 0) {
      vector.push_back(true);
      symbols.push_back(b);
    } else {
      vector.push_back(false);
    }
  }
  if(symbols.size() > 1) {
    outputShapeDfs(vector, 0, symbols);
  } else {
    vector.push_back(false);
  }
}

template <typename BitVector>
void WaveletTree<BitVector>::treeShape2(BitVector& vec) const {
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
  utils::pushBits(vec, packedInt, bytesInLongestCode*8);

  utils::binaryInterpolativeCode(symbols, symbols.back(), vec);
  for(size_t i = 0; i < symbols.size(); ++i) 
    utils::unaryCode(vec, maxLen - m_codes[symbols[i]].size() + 1);
}

template <typename BitVector>
void WaveletTree<BitVector>::outputShapeDfs(BitVector& output, size_t depth,
                                            const std::vector<byte>& symbols) const
{
  if (symbols.size() <= 1) return;
  std::vector<byte> leftSymbols, rightSymbols;
  for(size_t i = 0; i < symbols.size(); ++i) {
    if(m_codes[symbols[i]][depth]) {
      output.push_back(true);
      rightSymbols.push_back(symbols[i]);
    } else {
      output.push_back(false);
      leftSymbols.push_back(symbols[i]);
    }
  }
  outputShapeDfs(output, depth+1, leftSymbols);
  outputShapeDfs(output, depth+1, rightSymbols);
}

template <typename BitVector>
template <typename Encoder, typename ProbabilisticModel>
void WaveletTree<BitVector>::encodeTreeBF(Encoder& enc, ProbabilisticModel& pm) const
{
  std::queue<TreeNode<BitVector>*> queue;
  queue.push(m_root);
  while(!queue.empty()) {
    TreeNode<BitVector>* node = queue.front();
    for(size_t i = 0; i < node->m_bitVector.size(); ++i) {
      enc.encode(node->m_bitVector[i], pm.probabilityOfOne());
      pm.update(node->m_bitVector[i]);
    }
    queue.pop();
    if(node->m_left) queue.push(node->m_left);
    if(node->m_right) queue.push(node->m_right);
  }
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

/**Used in decoding tree.
 */
template<typename BitVector>
struct NodeStatus {
  //NodeStatus() : m_node(0), m_bits(0), m_gammaLen(0), m_gammaStatus(0) {}
  NodeStatus(TreeNode<BitVector>* node, size_t bits)
      : m_node(node), m_bits(bits), m_gammaLen(0), m_gammaStatus(0) {}
  NodeStatus(TreeNode<BitVector>* node, size_t bits, size_t gammaLen, byte st)
      : m_node(node), m_bits(bits), m_gammaLen(gammaLen), m_gammaStatus(st) {}
  TreeNode<BitVector> *m_node;
  size_t m_bits;
  size_t m_gammaLen;
  byte m_gammaStatus; // Flag telling if we have already hit into gamma codes
};

template<typename BitVector>
template <typename Decoder, typename ProbabilisticModel>
void WaveletTree<BitVector>::decodeTreeBF(size_t rootSize,
                                          Decoder& dec,
                                          ProbabilisticModel& pm) const
{
  std::queue<NodeStatus<BitVector> > queue;
  queue.push(NodeStatus<BitVector>(m_root, rootSize));
  while(!queue.empty()) {
    size_t rightBits = 0;
    NodeStatus<BitVector> st = queue.front();
    queue.pop();
    for(size_t i = 0; i < st.m_bits; ++i) {
      bool bit = dec.decode(pm.probabilityOfOne());
      pm.update(bit);
      if(bit) ++rightBits;
      st.m_node->m_bitVector.push_back(bit);
    }
    if(st.m_gammaStatus == 2 && st.m_gammaLen == 1) continue;

    if(st.m_bits > rightBits &&
       !(st.m_gammaStatus == 0 && st.m_node->m_hasSymbol))
    {
      if(!st.m_node->m_left) st.m_node->m_left = new TreeNode<BitVector>();
      NodeStatus<BitVector> left(st.m_node->m_left, st.m_bits - rightBits,
                                 st.m_gammaLen, st.m_gammaStatus);

      if(left.m_gammaStatus == 1) left.m_gammaStatus = 2;
      else if (left.m_gammaStatus == 2) --left.m_gammaLen;

      if(left.m_gammaStatus != 2 || left.m_gammaLen > 0) {
        queue.push(left);
      }
    }
    if(rightBits > 0) {
      if(!st.m_node->m_right) st.m_node->m_right = new TreeNode<BitVector>();
      NodeStatus<BitVector> right(st.m_node->m_right, rightBits,
                                  st.m_gammaLen, st.m_gammaStatus);

      if(st.m_gammaStatus == 0) 
        right.m_gammaStatus = st.m_node->m_hasSymbol ? 1 : 0;

      if(right.m_gammaStatus == 1) ++right.m_gammaLen;
      else if(right.m_gammaStatus == 2) --right.m_gammaLen;
        

      if(right.m_gammaStatus != 2 || right.m_gammaLen > 0) {
        queue.push(right);
      }

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

// TODO: optimize gamma coding
template <typename BitVector>
void WaveletTree<BitVector>::pushRun(byte symbol, size_t runLength)
{
  TreeNode<BitVector> *node = pushBits(m_codes[symbol]);
  assert(node->m_hasSymbol);
  BitVector integerCode;
  gammaCode(integerCode, runLength);
  pushBits(node, integerCode);
}

template <typename BitVector> template <typename OutputIterator>
void WaveletTree<BitVector>::message(OutputIterator out) const {
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
    // Decoding of gamma-code
    size_t runBits = 0;
    bit = node->m_bitVector[i];
    while(bit) {
      if(bit) i = bitsSeen[node]++;
      else i = i - bitsSeen[node]; 
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
    for(size_t k = 0; k < runLength; ++k) {
      *out++ = symbol;
    }
  }
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
    std::vector<std::pair<uint64, byte> >& lengths)
{
  std::sort(lengths.begin(), lengths.end());
  m_root = new TreeNode<BitVector>();
  assignPrefixCodes(lengths, m_root, 0, 0);
}

template <typename BitVector>
size_t WaveletTree<BitVector>::assignPrefixCodes(
    std::vector<std::pair<uint64, byte> >& lengths, TreeNode<BitVector>* node,
    size_t elem, size_t bits)
{
  if(elem >= lengths.size()) return elem;
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
void WaveletTree<BitVector>::collectCodes(BitVector *codes, TreeNode<BitVector> *root)
{
  BitVector vec;
  collectCodes(codes, vec, root);
}

template <typename BitVector>
void
WaveletTree<BitVector>::collectCodes(BitVector *codes, BitVector& vec, TreeNode<BitVector> *node)
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


} //namespace bwtc
#endif
