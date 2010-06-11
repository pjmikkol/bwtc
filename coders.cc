#include "coders.h"
#include "globaldefs.h"
#include "probmodels/base_prob_model.h"

#include <vector>

namespace bwtc {

ProbabilityModel* GiveProbabilityModel(char choice) {
  switch(choice) {
    case 'n':
    default:
        return new ProbabilityModel();
  }
}

Encoder::Encoder(const std::string& destination, char prob_model)
    : destination_(NULL), pm_(NULL) {
  destination_ = new dcsbwt::BitEncoder(destination);
  /* Add new probability models here: */
  pm_ = GiveProbabilityModel(prob_model);
}

void Encoder::EncodeByte(byte b) {
  for(int i = 0; i < 8; ++i, b <<= 1) {
    bool bit = b & 0x80;
    destination_->Encode(bit, pm_->ProbabilityOfOne());
    pm_->Update(bit);
  }
}

void Encoder::EncodeBlock(const byte* begin, const byte* end) {
  while(begin != end) {
    EncodeByte(*begin);
    ++begin;
  }
}

Encoder::~Encoder() {
  delete destination_;
  delete pm_;
}

Decoder::Decoder(const std::string& source, char prob_model) :
    source_(NULL), pm_(NULL) {
  source_ = new dcsbwt::BitDecoder(source);
  pm_ = GiveProbabilityModel(prob_model);
}

Decoder::~Decoder() {
  delete source_;
  delete pm_;
}

byte Decoder::DecodeByte() {
  byte b = 0x00;
  for(int i = 0; i < 8; ++i, b <<= 2) {
    if (source_->Decode(pm_->ProbabilityOfOne())) {
      b += 1;
      pm_->Update(true);
    } else {
      pm_->Update(false);
    }
  }
  return b;
}

} // namespace bwtc

