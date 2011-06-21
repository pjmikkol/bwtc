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

// Template argument is the number of different probabilities.
// The probabilites are conceptually in form ix, where i = 1,...,M-1 and x = 1/M.
template <uint16 M>
class BitPredictor : public ProbabilityModel {
 public:
  BitPredictor() : m_probabilityOfOne(m_interval * (M >> 2)) {}
  ~BitPredictor() {}

  void update(bool bit) {
    if(bit) {
      uint16 probability = m_probabilityOfOne + m_interval;
      if(probability < kProbabilityScale) m_probabilityOfOne = probability;
    } else {
      if(m_interval < m_probabilityOfOne) m_probabilityOfOne -= m_interval;
    }
  }

  Probability probabilityOfOne() const {
    return m_probabilityOfOne;
  }

  void resetModel() {
    m_probabilityOfOne = m_interval * (M >> 2);
  }
  
 private:
  static const uint16 m_interval = kProbabilityScale/M;
  uint16 m_probabilityOfOne;
  BOOST_STATIC_ASSERT(M <= kProbabilityScale);
  BOOST_STATIC_ASSERT(m_interval * M == kProbabilityScale);
  
};

} // namespace bwtc

#endif
