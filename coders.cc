#include <iostream> // For std::streampos
#include <numeric> // for std::accumulate
#include <vector>

#include "block.h"
#include "coders.h"
#include "globaldefs.h"
#include "probmodels/base_prob_model.h"

namespace bwtc {

uint64 PackInteger(uint64 integer, int* bytes_needed) {
  /* Results for the if-clause are undefined if integer value 0x80
   * doesn't have correct type (64-bit) */
  static const uint64 kEightBit = 0x80;

  uint64 result = 0; int i;
  // For optimization (if needed) two OR-operations could be merged
  for(i = 0; integer; ++i) {
    result |= ((integer & 0x7F) << i*8);
    integer >>= 7;
    assert(i < 8);
    if (integer) result |= (kEightBit << i*8);
  }
  *bytes_needed = i;
  return result;
}

uint64 UnpackInteger(uint64 packed_integer) {
  uint64 result = 0; int bits_handled = 0;
  bool bits_left;
  do {
    bits_left = (packed_integer & 0x80) != 0;
    result |= ((packed_integer & 0x7F) << bits_handled);
    packed_integer >>= 8;
    bits_handled += 7;
    assert(bits_handled <= 56);
  } while(bits_left);
  return result;
}


ProbabilityModel* GiveProbabilityModel(char choice) {
  switch(choice) {
    case 'n':
    default:
        return new ProbabilityModel();
  }
}

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
  assert(trailer);
  int bytes;
  uint64 packed_integer = PackInteger(trailer, &bytes);
  WritePackedInteger(packed_integer, bytes);
  return bytes;
}

/**************************************************************************
 *            Encoding and decoding single MainBlock                      *
 *------------------------------------------------------------------------*
 * Next functions handle encoding and decoding of the main blocks.        *
 * Also block header-format is specified here.                            *
 *                                                                        *
 * Format of the block is following:                                      *
 *  - Block header (no fixed length)                              (1)     *
 *  - Compressed block (no fixed length)                          (2)     *
 *  - Block trailer (coded same way as the context block lengths) (3)     *
 *                                                                        *
 * Block header format is following:                                      *
 * - Length of the header + compressed block + trailer in bytes (48 bits) *
 *   Note that the the length field itself isn't included                 *
 * - List of context block lengths in increasing order so that the stored *
 *   value is difference of current and previous length (always > 0).     *
 *   Lengths are compressed with PackInteger-function.                    *
 * - Ending symbol of the block header (2 bytes = 0x8000) which is        *
 *   invalid code for packed integer                                      *
 **************************************************************************/

// TODO: for compressing straight to the stdout we need to use
//       temporary file or huge buffer for each mainblock, so that
//       we can write the size of the compressed block in the beginning (1)
// At the moment the implementation is done only for compressing into file
void Encoder::EncodeData(std::vector<byte>* block, std::vector<uint64>* stats,
                         uint64 block_size) {
  // TODO: At the moment we are not assuming that block->frequencies_ is
  //       ordered in increasing order (since it isn't yet)
  unsigned i = 0;

  /* This loop is quite tricky since we want to prepare that context blocks
   * can be shattered around. */
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

void Encoder::WriteBlockHeader(std::vector<uint64>* stats) {
  uint64 header_length = 0;
  header_position_ = out_->GetPos();
  for (unsigned i = 0; i < 6; ++i) out_->WriteByte(0x00); //fill 48 bits
  for (unsigned i = 0; i < stats->size(); ++i) {
    int bytes;
    // TODO: At the moment we are not printing numbers in increasing order
    //       It has to be fixed at BWTransform and here
    if((*stats)[i] > 0) {
      uint64 packed_cblock_size = PackInteger((*stats)[i], &bytes);
      header_length += bytes;
      WritePackedInteger(packed_cblock_size, bytes);
    }
  }
  header_length += FinishBlockHeader();
  compressed_block_length_ = header_length;
  current_stat_handled_ = current_stat_index_ = 0;

  destination_->ResetCounter();
}

/* Integer is written in reversal fashion so that it can be read easier.*/
void Encoder::WritePackedInteger(uint64 packed_integer, int bytes) {
  do {
    byte to_written = static_cast<byte>(packed_integer & 0xFF);
    packed_integer >>= 8;
    out_->WriteByte(to_written);
  } while (--bytes);
  assert(0 == packed_integer);
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
    stats->push_back(UnpackInteger(value));
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
  *eof_byte = UnpackInteger(packed_integer);
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

