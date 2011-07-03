/**
 * @file FSMs.hpp
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
 * Templates for finite state automata-based probability models.
 */

#ifndef BWTC_FSMS_HPP
#define BWTC_FSMS_HPP

#include "../globaldefs.hpp"
#include "ProbabilityModel.hpp"

#include <algorithm>
#include <vector>
#include <boost/static_assert.hpp>

namespace bwtc {

namespace {

template<uint32 states>
uint32 nextState(uint32 currentState, bool bit) {
  if(bit) {
    if(currentState >= states/2)
      return std::min(currentState+1, states-1);
    else return states/2;
  } else {
    if(currentState < states/2)
      return std::max((int)currentState-1, 0);
    else return states/2 - 1;
  }
}

template<>
uint32 nextState<2>(uint32 currentState, bool bit) {
  return bit?1:0;
}

template<>
uint32 nextState<3>(uint32 currentState, bool bit) {
  if(currentState == 1) return bit?2:0;
  else if(currentState == 2 && bit) return 2;
  else if(currentState == 0 && !bit) return 0;
  else return 1;
}

}

/**Template parameters are: N for the number of states and BitPredictor for
 * predicting the bits inside state.
 */
template <uint32 N, typename BitPredictor>
class FSM : public ProbabilityModel {
 public:
  FSM() : m_currentState(N/2), m_states(N) {}
  ~FSM() {}

  void update(bool bit) {
    m_states[m_currentState].update(bit);
    m_currentState = nextState<N>(m_currentState, bit);
  }

  Probability probabilityOfOne() const {
    return m_states[m_currentState].probabilityOfOne();
  }

  void resetModel() {
    for(uint32 i = 0; i < N; ++i) m_states[i].resetModel();
    m_currentState = N/2;
  }

 private:
  uint32 m_currentState;
  std::vector<BitPredictor> m_states;
};

/** Each state has probability model attached:
 *  Z3 -- Z2 -- Z1 -- O1 -- O2 -- O3,
 *  where Z3 is used when at least 3 zeros are seen etc..
 */
template <typename Z3, typename Z2, typename Z1, typename O1, typename O2, typename O3>
class FSM6 : public ProbabilityModel {
 public:
  FSM6() : m_currentState(3) {}
  ~FSM6() {}

  void update(bool bit) {
    if(m_currentState <= 2) {
      if(m_currentState <= 1) {
        if (m_currentState == 0) z3.update(bit);
        else z2.update(bit);
      } else z1.update(bit);
    } else {
      if(m_currentState <= 4) {
        if (m_currentState == 3) o1.update(bit);
        else o2.update(bit);
      } else o3.update(bit);
    }
    m_currentState = nextState<6>(m_currentState, bit);
  }

  Probability probabilityOfOne() const {
    if(m_currentState <= 2) {
      if(m_currentState <= 1) {
        if (m_currentState == 0) return z3.probabilityOfOne();
        else return z2.probabilityOfOne();
      } else return z1.probabilityOfOne();
    } else {
      if(m_currentState <= 4) {
        if (m_currentState == 3) return o1.probabilityOfOne();
        else return o2.probabilityOfOne();
      } else return o3.probabilityOfOne();
    }
  }

  void resetModel() {
    z3.resetModel();
    z2.resetModel();
    z1.resetModel();
    o1.resetModel();
    o2.resetModel();
    o3.resetModel();
  }
  
  
 private:
  uint32 m_currentState;
  Z3 z3;
  Z2 z2;
  Z1 z1;
  O1 o1;
  O2 o2;
  O3 o3;
};

} // namespace bwtc

#endif
