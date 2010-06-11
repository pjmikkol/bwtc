#include "coders.h"
#include "probmodels/base_prob_model.h"

namespace bwtc {

Encoder::Encoder(const std::string& destination, char prob_model)
    : destination_(NULL), pm_(NULL) {
  destination_ = new dcsbwt::BitEncoder(destination);
  /* Add new probability models here: */
  pm_ = GiveProbabilityModel(prob_model);
}

Encoder::~Encoder() {
  delete destination_;
  delete pm_;
}

ProbabilityModel* GiveProbabilityModel(char choice) {
  switch(choice) {
    case 'n':
    default:
        return new ProbabilityModel();
  }
}

} // namespace bwtc

