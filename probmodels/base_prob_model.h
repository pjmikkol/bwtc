/* Base class acting as an interface for probability models. */
#ifndef BWTC_BASE_PROB_MODEL_H_
#define BWTC_BASE_PROB_MODEL_H_

#include "../globaldefs.h" /* Important definitions */

namespace bwtc {

class ProbabilityModel {
 public:
  ProbabilityModel() : prev_(true) {}
  virtual ~ProbabilityModel() {}
  virtual void Update(bool bit) { prev_ = bit; }
  virtual Probability ProbabilityOfOne() {
    if( prev_) return kProbabilityScale - 1;
    else return  1;
  }
  virtual void ResetModel() { prev_ = true; }

 private:
  bool prev_;
};

} // namespace bwtc

#endif
