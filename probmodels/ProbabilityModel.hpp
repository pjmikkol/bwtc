/**
 * @file ProbabilityModel.hpp
 * @author Pekka Mikkola <pjmikkol@cs.helsinki.fi>
 *
 * @section LICENSE
 *
 * This file is part of bwtc.
 *
 * bwtc is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * bwtc is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with bwtc.  If not, see <http://www.gnu.org/licenses/>.
 *
 * @section DESCRIPTION
 *
 * Header of the base class for probability models.
 */

/* Base class acting as an interface for probability models. */
#ifndef BWTC_PROBABILITY_MODEL_HPP_
#define BWTC_PROBABILITY_MODEL_HPP_

#include "../globaldefs.hpp" /* Important definitions */

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
  char* history_;
};

ProbabilityModel* GiveProbabilityModel(char choice);
  
} // namespace bwtc

#endif
