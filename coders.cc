/**************************************************************************
 *  Copyright 2010, Pekka Mikkola, pjmikkol (at) cs.helsinki.fi           *
 *                                                                        *
 *  This file is part of bwtc.                                            *
 *                                                                        *
 *  bwtc is free software: you can redistribute it and/or modify          *
 *  it under the terms of the GNU General Public License as published by  *
 *  the Free Software Foundation, either version 3 of the License, or     *
 *  (at your option) any later version.                                   *
 *                                                                        *
 *  bwtc is distributed in the hope that it will be useful,               *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *  GNU General Public License for more details.                          *
 *                                                                        *
 *  You should have received a copy of the GNU General Public License     *
 *  along with bwtc.  If not, see <http://www.gnu.org/licenses/>.         *
 **************************************************************************/

#include <cassert>

#include <iostream> // For std::streampos
#include <numeric> // for std::accumulate
#include <vector>

#include "block.h"
#include "coders.h"
#include "globaldefs.h"
#include "utils.h"
#include "probmodels/base_prob_model.h"

namespace bwtc {

Encoder::Encoder(const std::string& destination, char prob_model)
    : out_(NULL), destination_(NULL), pm_(NULL), header_position_(0),
      compressed_block_length_(0), current_stat_handled_(0),
      current_stat_index_(0)
{
  out_ = new OutStream(destination);
  destination_ = new dcsbwt::BitEncoder();
  destination_->Connect(out_);
  /* Add new probability models here: */
  pm_ = GiveProbabilityModel(prob_model);
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
  out_->WriteByte(static_cast<byte>(preproc));
  out_->WriteByte(static_cast<byte>(encoding));
}

char Decoder::ReadGlobalHeader() {
  char preproc = static_cast<char>(in_->ReadByte());
  char probmodel = static_cast<char>(in_->ReadByte());
  pm_ = GiveProbabilityModel(probmodel);
  return preproc;
}
/***************** Global header for file-section ends ******************/

void Encoder::EncodeByte(byte b) {
  for(int i = 0; i < 8; ++i, b <<= 1) {
    bool bit = b & 0x80;
    destination_->Encode(bit, pm_->ProbabilityOfOne());
    pm_->Update(bit);
  }
}

void Encoder::EncodeRange(const byte* begin, const byte* end) {
  while(begin != end) {
    EncodeByte(*begin);
    ++begin;
  }
}

Encoder::~Encoder() {
  delete out_;
  delete destination_;
  delete pm_;
}

void Encoder::EndContextBlock() {
  assert(pm_);
  pm_->ResetModel();
}

void Decoder::EndContextBlock() {
  assert(pm_);
  pm_->ResetModel();
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
  if (current_stat_handled_ == 0) 
    while((*stats)[current_stat_index_] == 0) ++current_stat_index_;
  
  for( ; ; ++current_stat_index_) {
    if (i == block_size) return;
    assert(current_stat_index_ < stats->size());
    for( ; current_stat_handled_ < (*stats)[current_stat_index_];
         ++i, ++current_stat_handled_) {
      if (i == block_size) return;
      EncodeByte((*block)[i]);
    }
    current_stat_handled_ = 0;
    EndContextBlock();
  }
}

void Encoder::FinishBlock(uint64 eob_byte) {
  destination_->Finish();
  compressed_block_length_ += destination_->Counter();
  compressed_block_length_ +=  WriteTrailer(eob_byte);
  out_->Write48bits(compressed_block_length_, header_position_);
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
  uint64 header_length = 0;
  header_position_ = out_->GetPos();
  for (unsigned i = 0; i < 6; ++i) out_->WriteByte(0x00); //fill 48 bits
  for (unsigned i = 0; i < stats->size(); ++i) {
    int bytes;
    // TODO: At the moment we are not printing numbers in increasing order
    //       It has to be fixed at BWTransform and here
    if((*stats)[i] > 0) {
      uint64 packed_cblock_size = utils::PackInteger((*stats)[i], &bytes);
      header_length += bytes;
      WritePackedInteger(packed_cblock_size);
    }
  }
  header_length += FinishBlockHeader();
  compressed_block_length_ = header_length;
  current_stat_handled_ = current_stat_index_ = 0;

  destination_->ResetCounter();
}

/* Integer is written in reversal fashion so that it can be read easier.*/
void Encoder::WritePackedInteger(uint64 packed_integer) {
  do {
    byte to_written = static_cast<byte>(packed_integer & 0xFF);
    packed_integer >>= 8;
    out_->WriteByte(to_written);
  } while (packed_integer);
}

int Encoder::FinishBlockHeader() {
  out_->WriteByte(0x80);
  out_->WriteByte(0x00);
  return 2;
}

uint64 Decoder::ReadBlockHeader(std::vector<uint64>* stats) {
  static const uint64 kErrorMask = static_cast<uint64>(1) << 63;

  uint64 compressed_length = in_->Read48bits();
  while(1) {
    uint64 value = ReadPackedInteger();
    if(value & kErrorMask) break;
    stats->push_back(utils::UnpackInteger(value));
  }
  return compressed_length;
}

std::vector<byte>* Decoder::DecodeBlock(uint64* eof_byte) {
  if(in_->CompressedDataEnding()) return NULL;

  std::vector<uint64> context_lengths;
  uint64 compr_len = ReadBlockHeader(&context_lengths);

  if (verbosity > 2)
    std::clog << "Size of compressed block = " << compr_len << "\n"; 

  uint64 block_size = std::accumulate(
      context_lengths.begin(), context_lengths.end(), static_cast<uint64>(0));
  source_->Start();
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
    uint64 read = static_cast<uint64>(in_->ReadByte());
    bits_left = (read & kEndMask) != 0;
    packed_integer |= (read << i*8);
  }
  if (packed_integer == 0x80) return kEndSymbol;
  return packed_integer;
}
/*********** Encoding and decoding single MainBlock-section ends ********/

Decoder::Decoder(const std::string& source, char prob_model) :
    in_(NULL), source_(NULL), pm_(NULL) {
  in_ = new InStream(source);
  source_ = new dcsbwt::BitDecoder();
  source_->Connect(in_);
  pm_ = GiveProbabilityModel(prob_model);
}

Decoder::Decoder(const std::string& source) :
    in_(NULL), source_(NULL), pm_(NULL) {
  in_ = new InStream(source);
  source_ = new dcsbwt::BitDecoder();
  source_->Connect(in_);
}

Decoder::~Decoder() {
  delete in_;
  delete source_;
  delete pm_;
}

byte Decoder::DecodeByte() {
  byte b = 0x00;
  for(int i = 0; i < 8; ++i) {
    b <<= 1;
    if (source_->Decode(pm_->ProbabilityOfOne())) {
      b |= 1;
      pm_->Update(true);
    } else {
      pm_->Update(false);
    }
  }
  return b;
}

} // namespace bwtc

