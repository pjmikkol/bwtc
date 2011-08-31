/**
 * @file BitPredictors.hpp
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
 * Templates for simple bit predictors.
 */

#ifndef BWTC_BIT_PREDICTORS_HPP
#define BWTC_BIT_PREDICTORS_HPP

#include "../globaldefs.hpp"
#include "ProbabilityModel.hpp"

#include <boost/static_assert.hpp>

namespace bwtc {

// TODO: LISÄÄ ALKUTILA
template <Probability Min, Probability Delay>
class UnbiasedPredictor : public ProbabilityModel {
 public:
  UnbiasedPredictor() { resetModel(); }
  ~UnbiasedPredictor() {}

  void update(bool bit) {
    if(bit) {
      m_probabilityOfOne += ((s_maximumProb - m_probabilityOfOne) >> Delay);
    } else {
      m_probabilityOfOne -= ((m_probabilityOfOne - Min) >> Delay);
    }
  }

  Probability probabilityOfOne() const {
    return m_probabilityOfOne;
  }
  
  void resetModel() {
    m_probabilityOfOne = kHalfProbability;
  }

 private:
  static const Probability s_maximumProb = kProbabilityScale - Min;
  Probability m_probabilityOfOne;

};


// Template argument is the number of different probabilities.
// The probabilites are conceptually in form ix, where i = 1,...,M-1 and x = 1/M.
template <Probability M>
class EvenIntervalPredictor : public ProbabilityModel {
 public:
  EvenIntervalPredictor() : m_probabilityOfOne(s_interval * (M >> 1)) {}
  ~EvenIntervalPredictor() {}

  void update(bool bit) {
    if(bit) {
      Probability probability = m_probabilityOfOne + s_interval;
      if(probability < kProbabilityScale) m_probabilityOfOne = probability;
    } else {
      if(s_interval < m_probabilityOfOne) m_probabilityOfOne -= s_interval;
    }
  }

  Probability probabilityOfOne() const {
    return m_probabilityOfOne;
  }

  void resetModel() {
    m_probabilityOfOne = s_interval * (M >> 1);
  }
  
 private:
  static const Probability s_interval = kProbabilityScale/M;
  Probability m_probabilityOfOne;
  BOOST_STATIC_ASSERT(M <= kProbabilityScale);
  BOOST_STATIC_ASSERT(s_interval * M == kProbabilityScale);
  
};

template <>
class EvenIntervalPredictor<2> : public ProbabilityModel {
 public:
  EvenIntervalPredictor() : m_probabilityOfOne(s_quarter + kHalfProbability) {}
  ~EvenIntervalPredictor() {}

  void update(bool bit) {
    if(bit) m_probabilityOfOne = s_quarter + kHalfProbability;
    else m_probabilityOfOne = s_quarter;
  }

  Probability probabilityOfOne() const {
    return m_probabilityOfOne;
  }

  void resetModel() {
    m_probabilityOfOne = s_quarter + kHalfProbability;
  }
  
 private:
  static const Probability s_quarter = kHalfProbability >> 1;
  Probability m_probabilityOfOne;
};

template <>
class EvenIntervalPredictor<3> : public ProbabilityModel {
 public:
  EvenIntervalPredictor() : m_probabilityOfOne(s_interval << 1) {}
  ~EvenIntervalPredictor() {}

  void update(bool bit) {
    if(bit) m_probabilityOfOne = s_interval << 1;
    else m_probabilityOfOne = s_interval;
  }

  Probability probabilityOfOne() const {
    return m_probabilityOfOne;
  }

  void resetModel() {
    m_probabilityOfOne = s_interval << 1;
  }
  
 private:
  static const Probability s_interval = kProbabilityScale/3;
  Probability m_probabilityOfOne;
};

template <>
class EvenIntervalPredictor<5> : public ProbabilityModel {
 public:
  EvenIntervalPredictor() : m_probabilityOfOne(s_interval << 1) {}
  ~EvenIntervalPredictor() {}

  void update(bool bit) {
    if(bit) {
      if (m_probabilityOfOne < s_interval)
        m_probabilityOfOne = s_interval;
      else if(m_probabilityOfOne < (s_interval << 1))
        m_probabilityOfOne = s_interval << 1;
      else
        m_probabilityOfOne = (s_interval << 1) + s_halfInterval;
    } else {
      if(m_probabilityOfOne > (s_interval << 1))
        m_probabilityOfOne = s_interval << 1;
      else if(m_probabilityOfOne > s_interval)
        m_probabilityOfOne = s_interval;
      else
        m_probabilityOfOne = s_interval - s_halfInterval;
    }
  }

  Probability probabilityOfOne() const {
    return m_probabilityOfOne;
  }

  void resetModel() {
    m_probabilityOfOne = s_interval << 1;
  }
  
 private:
  static const Probability s_interval = kProbabilityScale/3;
  static const Probability s_halfInterval = s_interval >> 1;
  Probability m_probabilityOfOne;
};


template <Probability M>
class BiasedOnePredictor : public ProbabilityModel {
 public:
  BiasedOnePredictor() : m_probabilityOfOne(s_interval * (M >> 1)) {}
  ~BiasedOnePredictor() {}

  void update(bool bit) {
    if(bit) {
      Probability increase = s_interval;
      Probability probability = m_probabilityOfOne + increase;

      while(probability >= kProbabilityScale) {
        probability -= increase;
        increase >>= 1;
        probability += increase;
      }
      m_probabilityOfOne = probability;
    } else {
      if((M-1)*s_interval < m_probabilityOfOne) m_probabilityOfOne = (M-1)*s_interval;
      else if(s_interval < m_probabilityOfOne) m_probabilityOfOne -= s_interval;
    }
  }

  Probability probabilityOfOne() const {
    return m_probabilityOfOne;
  }

  void resetModel() {
    m_probabilityOfOne = s_interval * (M >> 1);
  }
  
 private:
  static const Probability s_interval = kProbabilityScale/M;
  Probability m_probabilityOfOne;
  BOOST_STATIC_ASSERT(M <= kProbabilityScale);
  BOOST_STATIC_ASSERT(s_interval * M == kProbabilityScale);
  
};

class AggressiveOnePredictor : public ProbabilityModel {
 public:
  AggressiveOnePredictor() : m_probabilityOfOne(kHalfProbability), m_step(kHalfProbability) {}
  ~AggressiveOnePredictor() {}

  void update(bool bit) {
    if(bit) {
      if (m_probabilityOfOne < kHalfProbability) {
        m_probabilityOfOne = kHalfProbability;
        m_step = kHalfProbability;
      } else {
        if(m_step > 1) {
          m_step >>= 1;
          m_probabilityOfOne += m_step;
        }
      }
    } else {
      if(m_probabilityOfOne > kHalfProbability) {
        if (m_probabilityOfOne > s_limit) {
          m_probabilityOfOne = s_limit;
          m_step = kHalfProbability >> 3;
        } else {
          m_probabilityOfOne -= m_step;
          m_step <<= 1;
        }
      } else {
        m_probabilityOfOne = kHalfProbability >> 1;
      }
    }
  }

  Probability probabilityOfOne() const {
    return m_probabilityOfOne;
  }

  void resetModel() {
    m_probabilityOfOne = kHalfProbability;
  }

 private:
  static const Probability s_limit = kHalfProbability + (kHalfProbability >> 1) +
      (kHalfProbability >> 2) + (kHalfProbability >> 3);
  Probability m_probabilityOfOne;
  Probability m_step;
};

/** Inverses the predictions done by template predictor. */
template <typename Predictor>
class InversePredictor : public ProbabilityModel {
 public:
  InversePredictor() {}
  ~InversePredictor() {}

  void update(bool bit) {
    m_predictor.update(!bit);
  }

  Probability probabilityOfOne() const {
    return kProbabilityScale - m_predictor.probabilityOfOne();
  }

  void resetModel() {
    m_predictor.resetModel();
  }
  
 private:
  Predictor m_predictor;
};

} // namespace bwtc

#endif
