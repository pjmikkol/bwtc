/**
 * @file Coders.hpp
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
 * Header for coders.
 */

#ifndef BWTC_CODERS_HPP_
#define BWTC_CODERS_HPP_

#include <iostream>
#include <string>
#include <vector>

#include "MainBlock.hpp"
#include "globaldefs.hpp"
#include "BitCoders.hpp"
#include "probmodels/base_prob_model.h"

namespace bwtc {

/**********************************************************************
 * Encoder and decoder are pretty similar in structure.               *
 *                                                                    *
 * Both have a field for ProbabilityModel-object which ultimately     *
 * decides how the encoding or decoding is done.                      *
 *                                                                    *
 * Both have also destination or source field which is for arithmetic *
 * encoder/decoder (objects of a type BitEncoder/BitDecoder).         *
 **********************************************************************/

class Encoder {
 public:
  Encoder(const std::string& destination, char prob_model);
  ~Encoder();
  void WriteGlobalHeader(char preproc, char encoding);
  void EncodeByte(byte b);
  void EncodeData(std::vector<byte>* data, std::vector<uint64>* stats,
                  uint64 data_size);
  void EncodeRange(const byte* begin, const byte* end);
  void WriteBlockHeader(std::vector<uint64>* stats);
  void WritePackedInteger(uint64 packed_integer);
  int FinishBlockHeader();
  void EndContextBlock();
  int WriteTrailer(uint64 trailer_value);
  void FinishBlock(uint64 eob_byte); //TODO: this calls write trailer

 private:
  OutStream* m_out;
  dcsbwt::BitEncoder* m_destination;
  ProbabilityModel* m_probModel;
  std::streampos m_headerPosition;
  uint64 m_compressedBlockLength;
  /* We may have to encode result of the transformation in pieces so we
   * have track down the progress of handling single MainBlock. */
  uint64 m_currentStatHandled;
  unsigned m_currentStatIndex;


  Encoder(const Encoder&);
  Encoder& operator=(const Encoder&);
};

class Decoder {
 public:
  Decoder(const std::string& source, char prob_model);
  Decoder(const std::string& source);
  ~Decoder();
  /* ReadGlobalHeader returns char denoting the preprocessing algorithm.
   * It changes the used probability model automatically. */
  char ReadGlobalHeader();
  byte DecodeByte();
  void Start() { m_source->Start(); }
  /* If end symbol is encountered, then the most significant bit is activated */
  uint64 ReadPackedInteger();
  /* Allocates memory for block, reads and decodes it. */
  std::vector<byte>* DecodeBlock(uint64* eof_byte_in_bwt);
  /* Returns length of the compressed sequence and stores lengths of the context
   * blocks into stats-array.*/
  uint64 ReadBlockHeader(std::vector<uint64>* stats);
  void EndContextBlock();

 private:
  InStream* m_in;
  dcsbwt::BitDecoder* m_source;
  ProbabilityModel* m_probModel;

  Decoder(const Decoder&);
  Decoder& operator=(const Decoder&);
};

} // namespace bwtc

#endif
