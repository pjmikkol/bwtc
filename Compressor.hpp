/**
 * @file Compressor.hpp
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
 * Header for Compressor-class. The compressor is an abstraction of
 * an compression pipeline. It has 3 stages whic are illustrated in
 * the following diagram:
 *
 *  input -->  PRECOMPRESSION --> BWT --> ENTROPY CODING --> ouput
 *
 * The compressor has a memory limit which gives an upper bound for
 * the maximum amount of memory used. Limited by this restriction,
 * each phase handles as big blocks of data as possible.
 *
 * Precompression:
 * Precompressor reads as much data as possible into the memory and then
 * compresses the data read. After this precompressed data is divided into
 * BWT-blocks of maximum size (data has to be divided further since the
 * computation of BWT needs more memory than precompression).
 *
 * BWT:
 * Each BWT-block is transformed using BWTManager, which handles the
 * selection (and evaluation) of BWT-algorithms.
 *
 * Entropy coding:
 * Each BWT-block is compressed independently using some entropy coder.
 * Encoded BWTBlocks are written into the compressed file in the order
 * specified by ?TODO?
 *
 *
 * COMPRESSED FILE FORMAT:
 *
 * At the top level compressed file is divided into file header and several
 * precompression blocks.
 *
 * File header:
 * File header contains global information about the compressed file such as
 * used entropy coder etc.
 *
 * Precompression block:
 * Precompression blocks are independent of each other. Each precompression
 * block corresponds to a varying size of input data. Precompression block
 * contains header of the block and 1--n BWT-blocks. Header of Precompression
 * block contains metadata needed to uncompress precompression (in the case of
 * grammar compression that means grammar) and the number of BWT-blocks.
 *
 * BWT-block:
 * BWT-block contains header, trailer and entropy encoded data, which is
 * transformed. Header of BWT-block contains the size of the compressed BWT-
 * block and the size of uncompressed BWT-block (before precompression). In
 * addition header contains the data necessary to uncompress the block
 * (written by entropy encoder).
 * Trailer of BWT-block contains the number of starting points used in inverse
 * and their positions.
 *
 */

#ifndef BWTC_COMPRESSOR_HPP_
#define BWTC_COMPRESSOR_HPP_

#include "bwtransforms/BWTManager.hpp"
#include "preprocessors/Preprocessor.hpp"
#include "EntropyCoders.hpp"
#include "Streams.hpp"

#include <string>

namespace bwtc {

class Compressor {
 public:
  Compressor(const std::string& in, const std::string& out, size_t memLimit);
  Compressor(RawInStream* in, RawOutStream* out, size_t memLimit);
  ~Compressor();

  void setEntropyEncoder(char choice);
  void setPreprocessor(const std::string& parameters);

  void compress(size_t threads);

 private:
  RawInStream *m_in;
  RawOutStream *m_out;
  EntropyEncoder *m_coder;
  Preprocessor *m_preprocessor;
  //BWTManager m_bwtmanager;
  size_t m_memLimit;
};

} //namespace bwtc

#endif
