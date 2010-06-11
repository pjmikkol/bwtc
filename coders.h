#ifndef BWTC_CODERS_H_
#define BWTC_CODERS_H_

#include <string>

#include "globaldefs.h"
#include "rl_compress.h"
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

/* * 
 * In the case of encoding the ProbabilityModel should give probability
 * for each bit encoded. 
 * */

  //TODO: Decided to 
  //Should Encoder-class be a base class for different encoders
//      or should there be only a probability model field in Encoder-class
class Encoder {
 public:
  Encoder(const std::string& destination, char prob_model);
  ~Encoder();
  void EncodeByte(byte b);
  void EncodeBlock(const byte* begin, const byte* end);
  void Finish() { destination_->Finish(); }

 private:
  dcsbwt::BitEncoder* destination_;
  ProbabilityModel* pm_;

  Encoder(const Encoder&);
  Encoder& operator=(const Encoder&);
};

class Decoder {
 public:
  Decoder(const std::string& source, char prob_model);
  ~Decoder();
  byte DecodeByte();
  void Start() { source_->Start(); }

 private:
  dcsbwt::BitDecoder* source_;
  ProbabilityModel* pm_;

  Decoder(const Decoder&);
  Decoder& operator=(const Decoder&);
};

ProbabilityModel* GiveProbabilityModel(char choice);
  
} // namespace bwtc

#endif
