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

#include "EntropyCoders.hpp"
#include "MainBlock.hpp"
#include "globaldefs.hpp"
#include "BitCoders.hpp"
#include "probmodels/ProbabilityModel.hpp"
#include "Streams.hpp"
#include "BWTBlock.hpp"

namespace bwtc {

class WaveletEncoder : public EntropyEncoder {
 public:
  WaveletEncoder(char prob_model);
  ~WaveletEncoder();

  void encodeData(const byte* data, std::vector<uint64>* stats,
                  uint64 data_size, RawOutStream* out);
  void writeBlockHeader(std::vector<uint64>* stats, RawOutStream* out);
  void finishBlock(const std::vector<uint32>& LFpowers, RawOutStream* out);

  size_t transformAndEncode(BWTBlock& block, BWTManager& bwtm,
                            RawOutStream* out);
  
  void writePackedInteger(uint64 packed_integer, RawOutStream* out);
  void endContextBlock();
  int writeTrailer(const std::vector<uint32>& LFpowers, RawOutStream* out);

 private:
  dcsbwt::BitEncoder m_destination;
  /** Probability model for internal nodes in wavelet tree. */
  ProbabilityModel* m_probModel;
  /** Probability model for integer code nodes in wavelet tree. */
  ProbabilityModel* m_integerProbModel;
  /** Probability model for bits coming after gaps. */
  ProbabilityModel* m_gapProbModel;
  long int m_headerPosition;
  uint64 m_compressedBlockLength;
 
  WaveletEncoder(const WaveletEncoder&);
  WaveletEncoder& operator=(const WaveletEncoder&);
};

class WaveletDecoder : public EntropyDecoder {
 public:
  WaveletDecoder(const std::string& source);
  WaveletDecoder(RawInStream* in, char decoder);
  ~WaveletDecoder();
  /* If end symbol is encountered, then the most significant bit is activated */
  uint64 readPackedInteger();
  /* Allocates memory for block, reads and decodes it. */
  std::vector<byte>* decodeBlock(std::vector<uint32>& LFpowers);
  /* Returns length of the compressed sequence and stores lengths of the context
   * blocks into stats-array.*/
  uint64 readBlockHeader(std::vector<uint64>* stats);
  void endContextBlock();

 private:
  RawInStream* m_in;
  dcsbwt::BitDecoder m_source;
  /** Probability model for internal nodes in wavelet tree. */
  ProbabilityModel* m_probModel;
  /** Probability model for integer code nodes in wavelet tree. */
  ProbabilityModel* m_integerProbModel;
  /** Probability model for bits coming after gaps. */
  ProbabilityModel* m_gapProbModel;

  WaveletDecoder(const WaveletDecoder&);
  WaveletDecoder& operator=(const WaveletDecoder&);
};

} // namespace bwtc

#endif
