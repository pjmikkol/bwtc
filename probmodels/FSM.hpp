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

} // namespace bwtc

#endif
