/**
 * @file WaveletCoders.hpp
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
 * Header for coders which first construct wavelet tree and then compress it.
 */

#ifndef BWTC_WAVELET_CODERS_HPP_
#define BWTC_WAVELET_CODERS_HPP_

#include <iostream>
#include <string>
#include <vector>

#include "MainBlock.hpp"
#include "globaldefs.hpp"
#include "BitCoders.hpp"
#include "probmodels/ProbabilityModel.hpp"

namespace bwtc {

class WaveletEncoder {
 public:
  WaveletEncoder(const std::string& destination, char prob_model);
  ~WaveletEncoder();
  void writeGlobalHeader(char preproc, char encoding);
  void encodeData(std::vector<byte>* data, std::vector<uint64>* stats,
                  uint64 data_size);
  void writeBlockHeader(std::vector<uint64>* stats);

  void writePackedInteger(uint64 packed_integer);
  void endContextBlock();
  int writeTrailer(uint64 trailer_value);
  void finishBlock(uint64 eob_byte); //TODO: this calls write trailer

 private:
  OutStream* m_out;
  dcsbwt::BitEncoder m_destination;
  /** Probability model for internal nodes in wavelet tree. */
  ProbabilityModel* m_probModel;
  /** Probability model for gamma code nodes in wavelet tree. */
  ProbabilityModel* m_gammaProbModel;
  /** Probability model for bits coming after gaps. */
  ProbabilityModel* m_gapProbModel;
  std::streampos m_headerPosition;
  uint64 m_compressedBlockLength;

  WaveletEncoder(const WaveletEncoder&);
  WaveletEncoder& operator=(const WaveletEncoder&);
};

class WaveletDecoder {
 public:
  WaveletDecoder(const std::string& source);
  ~WaveletDecoder();
  /* ReadGlobalHeader returns char denoting the preprocessing algorithm.
   * It changes the used probability model automatically. */
  char readGlobalHeader();
  //void start() { m_source->start(); }
  /* If end symbol is encountered, then the most significant bit is activated */
  uint64 readPackedInteger();
  /* Allocates memory for block, reads and decodes it. */
  std::vector<byte>* decodeBlock(uint64* eof_byte_in_bwt);
  /* Returns length of the compressed sequence and stores lengths of the context
   * blocks into stats-array.*/
  uint64 readBlockHeader(std::vector<uint64>* stats);
  void endContextBlock();

 private:
  InStream* m_in;
  dcsbwt::BitDecoder m_source;
  /** Probability model for internal nodes in wavelet tree. */
  ProbabilityModel* m_probModel;
  /** Probability model for gamma code nodes in wavelet tree. */
  ProbabilityModel* m_gammaProbModel;
  /** Probability model for bits coming after gaps. */
  ProbabilityModel* m_gapProbModel;

  WaveletDecoder(const WaveletDecoder&);
  WaveletDecoder& operator=(const WaveletDecoder&);
};

} // namespace bwtc

#endif
