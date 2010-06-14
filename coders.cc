#include <iostream> // For std::streamsize

#include "block.h"
#include "coders.h"
#include "globaldefs.h"
#include "probmodels/base_prob_model.h"

namespace bwtc {

int64 PackInteger(int64 integer, int* bytes_needed) {
  return 1;
}

int64 UnpackInteger(int64 packed_integer) {
  return 1;
}


ProbabilityModel* GiveProbabilityModel(char choice) {
  switch(choice) {
    case 'n':
    default:
        return new ProbabilityModel();
  }
}

Encoder::Encoder(const std::string& destination, char prob_model)
    : out_(NULL), destination_(NULL), pm_(NULL) {
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

/*************************************************************************
 *            Encoding and decoding single MainBlock                     *
 *-----------------------------------------------------------------------*
 * Next functions handle encoding and decoding of the main blocks.       *
 * Also block header-format is specified here.                           *
 *************************************************************************/

// TODO: for compressing straight to the stdout we need to use
//       temporary file or huge buffer for each mainblock, so that
//       we can write the size of the compressed block in the beginning
// At the moment the implementation is done only for compressing into file
void Encoder::EncodeMainBlock(bwtc::MainBlock* block) {
  std::streampos len_pos = WriteBlockHeader(block->Stats());
}

std::streampos Encoder::WriteBlockHeader(int64* stats) {
  std::streampos header_start = out_->GetPos();
  for (int i = 0; i < 6; ++i) out_->WriteByte(0x00); //fill 48 bits
  for (int i = 0; i < 256; ++i) {
    int bytes;
    if(stats[i]) {
      int64 packed_cblock_size = PackInteger(stats[i], &bytes);
    }
  }
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
    if (source_->Decode(pm_->ProbabilityOfOne())) {
      b |= 1;
      pm_->Update(true);
    } else {
      pm_->Update(false);
    }
    if (i < 7) b <<= 1;
  }
  return b;
}

} // namespace bwtc

