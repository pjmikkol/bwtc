/**
 * @file Coders.cpp
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
 * Implementations for the decoder and encoder.
 */

#include <cassert>

#include <iostream> // For std::streampos
#include <numeric> // for std::accumulate
#include <vector>

#include "MainBlock.hpp"
#include "Coders.hpp"
#include "globaldefs.hpp"
#include "Utils.hpp"
#include "probmodels/base_prob_model.h"

namespace bwtc {

Encoder::Encoder(const std::string& destination, char prob_model)
    : m_out(new OutStream(destination)), m_destination(new dcsbwt::BitEncoder()),
      m_probModel(GiveProbabilityModel(prob_model)), m_headerPosition(0),
      m_compressedBlockLength(0), m_currentStatHandled(0),
      m_currentStatIndex(0)
{
  m_destination->Connect(m_out);
}

/************************************************************************
 *                Global header for file                                *
 *----------------------------------------------------------------------*
 * Write- and ReadGlobalHeader defines global header and its format for *
 * the compressed file.                                                 *
 ************************************************************************/
void Encoder::WriteGlobalHeader(char preproc, char encoding) {
  /* At the moment dummy implementation. In future should use
   * bit-fields of a bytes as a flags. */
  m_out->WriteByte(static_cast<byte>(preproc));
  m_out->WriteByte(static_cast<byte>(encoding));
}

char Decoder::ReadGlobalHeader() {
  char preproc = static_cast<char>(m_in->ReadByte());
  char probmodel = static_cast<char>(m_in->ReadByte());
  m_probModel = GiveProbabilityModel(probmodel);
  return preproc;
}
/***************** Global header for file-section ends ******************/

void Encoder::EncodeByte(byte b) {
  for(int i = 0; i < 8; ++i, b <<= 1) {
    bool bit = b & 0x80;
    m_destination->Encode(bit, m_probModel->ProbabilityOfOne());
    m_probModel->Update(bit);
  }
}

void Encoder::EncodeRange(const byte* begin, const byte* end) {
  while(begin != end) {
    EncodeByte(*begin);
    ++begin;
  }
}

Encoder::~Encoder() {
  delete m_out;
  delete m_destination;
  delete m_probModel;
}

void Encoder::EndContextBlock() {
  assert(m_probModel);
  m_probModel->ResetModel();
}

void Decoder::EndContextBlock() {
  assert(m_probModel);
  m_probModel->ResetModel();
}

int Encoder::WriteTrailer(uint64 trailer) {
  int bytes;
  uint64 packed_integer = utils::PackInteger(trailer, &bytes);
  WritePackedInteger(packed_integer);
  return bytes;
}

/******************************************************************************
 *            Encoding and decoding single MainBlock                          *
 *----------------------------------------------------------------------------*
 * Following functions handle encoding and decoding of the main blocks.       *
 * Also block header-format is specified here.                                *
 *                                                                            *
 * Format of the block is following:                                          *
 *  - Block header (no fixed length)                              (1)         *
 *  - Compressed block (no fixed length)                          (2)         *
 *  - Block trailer (coded same way as the context block lengths) (3)         *
 *                                                                            *
 *----------------------------------------------------------------------------*
 *                                                                            *
 * Block header (1) format is following:                                      *
 * a) Length of the (header + compressed block + trailer) in bytes (6 bytes). *
 *    Note that the the length field itself isn't includedin total length     *
 * b) List of context block lengths. Lengths are compressed with              *
 *    utils::PackInteger-function.                                            *
 * c) Ending symbol of the block header (2 bytes = 0x8000) which is           *
 *    invalid code for packed integer                                         *
 ******************************************************************************/

// TODO: for compressing straight to the stdout we need to use
//       temporary file or huge buffer for each mainblock, so that
//       we can write the size of the compressed block in the beginning (1a)
// At the moment the implementation is done only for compressing into file
void Encoder::EncodeData(std::vector<byte>* block, std::vector<uint64>* stats,
                         uint64 block_size)
{
  /* At the moment the context block-lengths are written in same order than   *
   * are found form the stats. Transform which build stats, have to re-order  *
   * lengths of stat-array if we are going to use some more clever packing of *
   * the context block-lengths (eg. ascending/descending order)               */
  unsigned i = 0;
  /* This loop is quite tricky since we want to be prepared that context
   * blocks can be shattered around. */
  if (m_currentStatHandled == 0) 
    while((*stats)[m_currentStatIndex] == 0) ++m_currentStatIndex;
  
  for( ; ; ++m_currentStatIndex) {
    if (i == block_size) return;
    assert(m_currentStatIndex < stats->size());
    for( ; m_currentStatHandled < (*stats)[m_currentStatIndex];
         ++i, ++m_currentStatHandled) {
      if (i == block_size) return;
      EncodeByte((*block)[i]);
    }
    m_currentStatHandled = 0;
    EndContextBlock();
  }
}

void Encoder::FinishBlock(uint64 eob_byte) {
  m_destination->Finish();
  m_compressedBlockLength += m_destination->Counter();
  m_compressedBlockLength +=  WriteTrailer(eob_byte);
  m_out->Write48bits(m_compressedBlockLength, m_headerPosition);
}

/*********************************************************************
 * The format of header for single main block is the following:      *
 * - 48 bits for the length of the compressed main block, doesn't    *
 *   include 6 bytes used for this                                   *
 * - lengths of context blocks in the same order as they are in the  *
 *   result of BWT. Lengths are coded with utils::PackInteger        *
 * - 2 sentinel bytes 0x80 and 0x00 for notifying end of the header  *
 *********************************************************************/
void Encoder::WriteBlockHeader(std::vector<uint64>* stats) {
  uint64 headerLength = 0;
  m_headerPosition = m_out->GetPos();
  for (unsigned i = 0; i < 6; ++i) m_out->WriteByte(0x00); //fill 48 bits
  for (unsigned i = 0; i < stats->size(); ++i) {
    int bytes;
    // TODO: At the moment we are not printing numbers in increasing order
    //       It has to be fixed at BWTransform and here
    if((*stats)[i] > 0) {
      uint64 packed_cblock_size = utils::PackInteger((*stats)[i], &bytes);
      headerLength += bytes;
      WritePackedInteger(packed_cblock_size);
    }
  }
  headerLength += FinishBlockHeader();
  m_compressedBlockLength = headerLength;
  m_currentStatHandled = m_currentStatIndex = 0;

  m_destination->ResetCounter();
}

/* Integer is written in reversal fashion so that it can be read easier.*/
void Encoder::WritePackedInteger(uint64 packed_integer) {
  do {
    byte to_written = static_cast<byte>(packed_integer & 0xFF);
    packed_integer >>= 8;
    m_out->WriteByte(to_written);
  } while (packed_integer);
}

int Encoder::FinishBlockHeader() {
  m_out->WriteByte(0x80);
  m_out->WriteByte(0x00);
  return 2;
}

uint64 Decoder::ReadBlockHeader(std::vector<uint64>* stats) {
  static const uint64 kErrorMask = static_cast<uint64>(1) << 63;

  uint64 compressed_length = m_in->Read48bits();
  while(1) {
    uint64 value = ReadPackedInteger();
    if(value & kErrorMask) break;
    stats->push_back(utils::UnpackInteger(value));
  }
  return compressed_length;
}

std::vector<byte>* Decoder::DecodeBlock(uint64* eof_byte) {
  if(m_in->CompressedDataEnding()) return NULL;

  std::vector<uint64> context_lengths;
  uint64 compr_len = ReadBlockHeader(&context_lengths);

  if (verbosity > 2) {
    std::clog << "Size of compressed block = " << compr_len << "\n";
  }

  uint64 block_size = std::accumulate(
      context_lengths.begin(), context_lengths.end(), static_cast<uint64>(0));
  m_source->Start();
  std::vector<byte>* data = new std::vector<byte>(block_size);
  int j = 0;
  for(std::vector<uint64>::const_iterator it = context_lengths.begin();
      it != context_lengths.end(); ++it) {
    for(uint64 i = 0; i < *it; ++i) {
      (*data)[j++] = DecodeByte();
    }
    EndContextBlock();
  }
  uint64 packed_integer = ReadPackedInteger();
  *eof_byte = utils::UnpackInteger(packed_integer);
  return data;
}

uint64 Decoder::ReadPackedInteger() {
  static const uint64 kEndSymbol = static_cast<uint64>(1) << 63;
  static const uint64 kEndMask = static_cast<uint64>(1) << 7;

  uint64 packed_integer = 0;
  bool bits_left = true;
  int i;
  for(i = 0; bits_left; ++i) {
    uint64 read = static_cast<uint64>(m_in->ReadByte());
    bits_left = (read & kEndMask) != 0;
    packed_integer |= (read << i*8);
  }
  if (packed_integer == 0x80) return kEndSymbol;
  return packed_integer;
}
/*********** Encoding and decoding single MainBlock-section ends ********/

Decoder::Decoder(const std::string& source, char prob_model)
    : m_in(new InStream(source)), m_source(new dcsbwt::BitDecoder()),
      m_probModel(GiveProbabilityModel(prob_model))
{
  m_source->Connect(m_in);
}

Decoder::Decoder(const std::string& source) :
    m_in(new InStream(source)), m_source(new dcsbwt::BitDecoder()),
    m_probModel(0)
{
  m_source->Connect(m_in);
}

Decoder::~Decoder() {
  delete m_in;
  delete m_source;
  delete m_probModel;
}

byte Decoder::DecodeByte() {
  byte b = 0x00;
  for(int i = 0; i < 8; ++i) {
    b <<= 1;
    if (m_source->Decode(m_probModel->ProbabilityOfOne())) {
      b |= 1;
      m_probModel->Update(true);
    } else {
      m_probModel->Update(false);
    }
  }
  return b;
}

} // namespace bwtc

