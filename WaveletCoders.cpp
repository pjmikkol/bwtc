/**
 * @file WaveletCoders.cpp
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
 * Implementations for the wavelet-decoder and -encoder.
 */

#include <cassert>

#include <iostream> // For std::streampos
#include <numeric> // for std::accumulate
#include <vector>

#include "MainBlock.hpp"
#include "Coders.hpp"
#include "globaldefs.hpp"
#include "Utils.hpp"
#include "probmodels/ProbabilityModel.hpp"

namespace bwtc {

WaveletEncoder::WaveletEncoder(const std::string& destination, char prob_model)
    : m_out(new OutStream(destination)),
      m_probModel(GiveProbabilityModel(prob_model)),
      m_headerPosition(0), m_compressedBlockLength(0), m_currentStatHandled(0),
      m_currentStatIndex(0)
{
  m_destination.connect(m_out);
}

void WaveletEncoder::writeGlobalHeader(char preproc, char encoding) {
  /* At the moment dummy implementation. In future should use
   * bit-fields of a bytes as a flags. */
  m_out->writeByte(static_cast<byte>(preproc));
  m_out->writeByte(static_cast<byte>(encoding));
}

char WaveletDecoder::readGlobalHeader() {
  char preproc = static_cast<char>(m_in->readByte());
  char probmodel = static_cast<char>(m_in->readByte());
  delete m_probModel;
  m_probModel = giveProbabilityModel(probmodel);
  return preproc;
}

WaveletEncoder::~WaveletEncoder() {
  delete m_out;
  delete m_probModel;
}

void WaveletEncoder::endContextBlock() {
  assert(m_probModel);
  m_probModel->resetModel();
}

void WaveletDecoder::endContextBlock() {
  assert(m_probModel);
  m_probModel->resetModel();
}

int WaveletEncoder::writeTrailer(uint64 trailer) {
  int bytes;
  uint64 packed_integer = utils::packInteger(trailer, &bytes);
  writePackedInteger(packed_integer);
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
 *    Note that the the length field itself isn't included in total length    *
 * b) List of context block lengths. Lengths are compressed with              *
 *    utils::PackInteger-function.                                            *
 * c) Ending symbol of the block header (2 bytes = 0x8000) which is           *
 *    invalid code for packed integer                                         *
 ******************************************************************************/

// TODO: for compressing straight to the stdout we need to use
//       temporary file or huge buffer for each mainblock, so that
//       we can write the size of the compressed block in the beginning (1a)
// At the moment the implementation is done only for compressing into file
void WaveletEncoder::encodeData(std::vector<byte>* block, std::vector<uint64>* stats,
                                uint64 block_size)
{
  size_t beg = 0;
  for(size_t i = 0; i < stats->size(); ++i) {
    if((*stats)[i] == 0) continue;
    WaveletTree wavelet(&(*block)[beg], (*stats)[i] + beg);
    // TODO: compression of tree
    beg += (*stats)[i];
    endContextBlock();
  }
}

void WaveletEncoder::finishBlock(uint64 eob_byte) {
  m_destination.finish();
  m_compressedBlockLength += m_destination.counter();
  m_compressedBlockLength +=  writeTrailer(eob_byte);
  m_out->write48bits(m_compressedBlockLength, m_headerPosition);
}

/*********************************************************************
 * The format of header for single main block is the following:      *
 * - 48 bits for the length of the compressed main block, doesn't    *
 *   include 6 bytes used for this                                   *
 * - lengths of context blocks in the same order as they are in the  *
 *   result of BWT. Lengths are coded with utils::PackInteger        *
 * - 2 sentinel bytes 0x80 and 0x00 for notifying end of the header  *
 *********************************************************************/
void WaveletEncoder::writeBlockHeader(std::vector<uint64>* stats) {
  uint64 headerLength = 0;
  m_headerPosition = m_out->getPos();
  for (unsigned i = 0; i < 6; ++i) m_out->writeByte(0x00); //fill 48 bits
  for (unsigned i = 0; i < stats->size(); ++i) {
    int bytes;
    // TODO: Maybe more clever encoding
    if((*stats)[i] > 0) {
      uint64 packed_cblock_size = utils::packInteger((*stats)[i], &bytes);
      headerLength += bytes;
      writePackedInteger(packed_cblock_size);
    }
  }
  headerLength += finishBlockHeader();
  m_compressedBlockLength = headerLength;
  m_currentStatHandled = m_currentStatIndex = 0;

  m_destination.resetCounter();
}

/* Integer is written in reversal fashion so that it can be read easier.*/
void WaveletEncoder::writePackedInteger(uint64 packed_integer) {
  do {
    byte to_written = static_cast<byte>(packed_integer & 0xFF);
    packed_integer >>= 8;
    m_out->writeByte(to_written);
  } while (packed_integer);
}

int WaveletEncoder::finishBlockHeader() {
  m_out->writeByte(0x80);
  m_out->writeByte(0x00);
  return 2;
}

uint64 WaveletDecoder::readBlockHeader(std::vector<uint64>* stats) {
  static const uint64 kErrorMask = static_cast<uint64>(1) << 63;

  uint64 compressed_length = m_in->read48bits();
  while(1) {
    uint64 value = readPackedInteger();
    if(value & kErrorMask) break;
    stats->push_back(utils::unpackInteger(value));
  }
  return compressed_length;
}

std::vector<byte>* WaveletDecoder::decodeBlock(uint64* eof_byte) {
  if(m_in->compressedDataEnding()) return NULL;

  std::vector<uint64> context_lengths;
  uint64 compr_len = readBlockHeader(&context_lengths);

  if (verbosity > 2) {
    std::clog << "Size of compressed block = " << compr_len << "\n";
  }

  uint64 block_size = std::accumulate(
      context_lengths.begin(), context_lengths.end(), static_cast<uint64>(0));
  m_source.start();
  std::vector<byte>* data = new std::vector<byte>(block_size);
  int j = 0;
  for(std::vector<uint64>::const_iterator it = context_lengths.begin();
      it != context_lengths.end(); ++it) {
    for(uint64 i = 0; i < *it; ++i) {
      //TODO: uncompress wavelet tree
    }
    endContextBlock();
  }
  uint64 packed_integer = readPackedInteger();
  *eof_byte = utils::unpackInteger(packed_integer);
  return data;
}

uint64 WaveletDecoder::readPackedInteger() {
  static const uint64 kEndSymbol = static_cast<uint64>(1) << 63;
  static const uint64 kEndMask = static_cast<uint64>(1) << 7;

  uint64 packed_integer = 0;
  bool bits_left = true;
  int i;
  for(i = 0; bits_left; ++i) {
    uint64 read = static_cast<uint64>(m_in->readByte());
    bits_left = (read & kEndMask) != 0;
    packed_integer |= (read << i*8);
  }
  if (packed_integer == 0x80) return kEndSymbol;
  return packed_integer;
}
/*********** Encoding and decoding single MainBlock-section ends ********/

WaveletDecoder::WaveletDecoder(const std::string& source, char prob_model)
    : m_in(new InStream(source)),
      m_probModel(GiveProbabilityModel(prob_model))
{
  m_source.connect(m_in);
}

WaveletDecoder::WaveletDecoder(const std::string& source) :
    m_in(new InStream(source)), m_probModel(0)
{
  m_source.connect(m_in);
}

WaveletDecoder::~WaveletDecoder() {
  delete m_in;
  delete m_probModel;
}


} // namespace bwtc

