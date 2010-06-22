/* Base class acting as an interface for probability models. */
#ifndef BWTC_BASE_PROB_MODEL_H_
#define BWTC_BASE_PROB_MODEL_H_

#include "../globaldefs.h" /* Important definitions */

namespace bwtc {

class ProbabilityModel {
 public:
  ProbabilityModel() : prev_(true) {}
  virtual ~ProbabilityModel() {}
  /* This will be called each time after single bit is coded. Updates to
   * model should be done here. */
  virtual void Update(bool bit) { prev_ = bit; }
  /* This probability will be used for coding each bit of the source. */
  virtual Probability ProbabilityOfOne() {
    if( prev_) return kProbabilityScale - 1;
    else return  1;
  }
  /* This will called when the context of data changes. */
  virtual void ResetModel() { prev_ = true; }

 private:
  bool prev_;
};

/* Example how to integrate new probability model to program */
template <typename UnsignedInt>
class SimpleMarkov : public ProbabilityModel {
 public:
  SimpleMarkov();
  virtual ~SimpleMarkov();
  virtual void Update(bool bit);
  virtual Probability ProbabilityOfOne();
  virtual void ResetModel();

 private:
  UnsignedInt prev_;
  bool* history_;
};

ProbabilityModel* GiveProbabilityModel(char choice);
  
} // namespace bwtc

#endif
