#include "encoder.h"

namespace bwtc {

Encoder* GiveEncoder(char choice, const std::string& destination) {
  /* Expand this to conform for different encoding algorithms */
  Encoder* enc;
  switch(choice) {
    case 'n':
    default:
        enc = new Encoder();
  }
  enc->Connect(destination);
  return enc;
}

Encoder::Encoder() : destination_(NULL) {}

Encoder::~Encoder() {
  delete destination_;
}

void Encoder::Connect(const std::string& destination) {
  destination_ = new dcsbwt::BitEncoder(destination);
}

} // namespace bwtc

