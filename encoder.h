#include <string>

#include "rl_compress.h"

namespace bwtc {

//TODO: Should Encoder-class be a base class for different encoders
//      or should there be only a probability model field in Encoder-class


class Encoder {
 public:
  Encoder();
  ~Encoder();
  void Connect(const std::string& destination); 

 private:
  dcsbwt::BitEncoder* destination_;

  Encoder(const Encoder&);
  Encoder& operator=(const Encoder&);
};

Encoder* GiveEncoder(char choice, const std::string& destination);


} // namespace bwtc
