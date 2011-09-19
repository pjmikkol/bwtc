/**
 * @file FSM.hpp
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
#include "BitPredictors.hpp"

#include <algorithm>
#include <vector>

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
    else return (states-1)/2;
  }
}

template<>
uint32 nextState<2>(uint32 , bool bit) {
  return bit?1:0;
}

template<>
uint32 nextState<3>(uint32 currentState, bool bit) {
  if(currentState == 1) return bit?2:0;
  else if(currentState == 2 && bit) return 2;
  else if(currentState == 0 && !bit) return 0;
  else return 1;
}

template<>
uint32 nextState<9>(uint32 currentState, bool bit) {
  static byte oneTransitions[] = {5, 5, 5, 4, 5, 6, 7, 8, 8};
  if(bit) return oneTransitions[currentState];
  else return 8 - oneTransitions[8 - currentState];
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
    updateState(bit);
  }

  Probability probabilityOfOne() const {
    return m_states[m_currentState].probabilityOfOne();
  }

  void resetModel() {
    for(uint32 i = 0; i < N; ++i) m_states[i].resetModel();
    m_currentState = N/2;
  }

  void updateState(bool bit) {
    m_currentState = nextState<N>(m_currentState, bit); 
  }

 private:
  uint32 m_currentState;
  std::vector<BitPredictor> m_states;
};

/** Each state has probability model attached:
 *  z3 -- z2 -- z1 -- o1 -- o2 -- o3,
 *  where Z3 is used when at least 3 zeros are seen etc..
 *  Types of the models are Z3, Z2, Z1, Z1, Z2, Z3.
 */
template <typename Z3, typename Z2, typename Z1>
class FSM6 : public ProbabilityModel {
 public:
  FSM6() : m_currentState(3) {}
  ~FSM6() {}

  void update(bool bit) {
    switch(m_currentState) {
      case 0: z3.update(bit); break;
      case 1: z2.update(bit); break;
      case 2: z1.update(bit); break;
      case 3: o1.update(bit); break;
      case 4: o2.update(bit); break;
      case 5: o3.update(bit); break;
    }
    updateState(bit);
  }

  void updateState(bool bit) {
    m_currentState = nextState<6>(m_currentState, bit);
  }

  Probability probabilityOfOne() const {
    switch(m_currentState) {
      case 0: return z3.probabilityOfOne();
      case 1: return z2.probabilityOfOne();
      case 2: return z1.probabilityOfOne();
      case 3: return o1.probabilityOfOne();
      case 4: return o2.probabilityOfOne();
      default: case 5: return o3.probabilityOfOne();
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
  InversePredictor<Z1> o1;
  InversePredictor<Z2> o2;
  InversePredictor<Z3> o3;
};

/** Almost the same as FSM6.
 */
template <typename Z4, typename Z3, typename Z2, typename Z1>
class FSM8 : public ProbabilityModel {
 public:
  FSM8() : m_currentState(4) {}
  ~FSM8() {}

  void update(bool bit) {
    switch(m_currentState) {
      case 0: z4.update(bit); break;
      case 1: z3.update(bit); break;
      case 2: z2.update(bit); break;
      case 3: z1.update(bit); break;
      case 4: o1.update(bit); break;
      case 5: o2.update(bit); break;
      case 6: o3.update(bit); break;
      case 7: o4.update(bit); break;
    }
    updateState(bit);
  }

  void updateState(bool bit) {
    m_currentState = nextState<8>(m_currentState, bit);
  }    
  
  Probability probabilityOfOne() const {
    switch(m_currentState) {
      case 0: return z4.probabilityOfOne();
      case 1: return z3.probabilityOfOne();
      case 2: return z2.probabilityOfOne();
      case 3: return z1.probabilityOfOne();
      case 4: return o1.probabilityOfOne();
      case 5: return o2.probabilityOfOne();
      case 6: return o3.probabilityOfOne();
      default: case 7: return o4.probabilityOfOne();
    }
  }

  void resetModel() {
    z4.resetModel();
    z3.resetModel();
    z2.resetModel();
    z1.resetModel();
    o1.resetModel();
    o2.resetModel();
    o3.resetModel();
    o4.resetModel();
  }
  
  
 private:
  uint32 m_currentState;
  Z4 z4;
  Z3 z3;
  Z2 z2;
  Z1 z1;
  InversePredictor<Z1> o1;
  InversePredictor<Z2> o2;
  InversePredictor<Z3> o3;
  InversePredictor<Z4> o4;
};

template <typename Z4, typename Z3, typename Z2, typename Z1, typename D>
class FSM9 : public ProbabilityModel {
 public:
  FSM9() : m_currentState(4) {}
  ~FSM9() {}

  void update(bool bit) {
    switch(m_currentState) {
      case 0: z4.update(bit); break;
      case 1: z3.update(bit); break;
      case 2: z2.update(bit); break;
      case 3: z1.update(bit); break;
      case 4: d.update(bit); break;
      case 5: o1.update(bit); break;
      case 6: o2.update(bit); break;
      case 7: o3.update(bit); break;
      case 8: o4.update(bit); break;
    }
    updateState(bit);
  }

  void updateState(bool bit) {
    m_currentState = nextState<9>(m_currentState, bit);
  }

  Probability probabilityOfOne() const {
    switch(m_currentState) {
      case 0: return z4.probabilityOfOne();
      case 1: return z3.probabilityOfOne();
      case 2: return z2.probabilityOfOne();
      case 3: return z1.probabilityOfOne();
      case 4: return d.probabilityOfOne();
      case 5: return o1.probabilityOfOne();
      case 6: return o2.probabilityOfOne();
      case 7: return o3.probabilityOfOne();
      default: case 8: return o4.probabilityOfOne();
    }
  }

  void resetModel() {
    z4.resetModel();
    z3.resetModel();
    z2.resetModel();
    z1.resetModel();
    d.resetModel();
    o1.resetModel();
    o2.resetModel();
    o3.resetModel();
    o4.resetModel();
  }

 private:
  uint32 m_currentState;
  Z4 z4;
  Z3 z3;
  Z2 z2;
  Z1 z1;
  D d;
  InversePredictor<Z1> o1;
  InversePredictor<Z2> o2;
  InversePredictor<Z3> o3;
  InversePredictor<Z4> o4;
};

template <typename BitPredictor, uint32 History>
class LimitedHistoryModel : public ProbabilityModel {
 public:
  LimitedHistoryModel() : m_history(0xAAAAAAAA & s_mask) {
    m_predictors.resize(1 << History);
  }

  ~LimitedHistoryModel() {}

  void update(bool bit) {
    m_predictors[m_history].update(bit);
    updateState(bit);
  }
  
  Probability probabilityOfOne() const {
    return m_predictors[m_history].probabilityOfOne();
  }
  
  void updateState(bool bit) {
    m_history = ((m_history << 1) + (bit?1:0))&s_mask;
  }
  

 private:
  std::vector<BitPredictor> m_predictors;
  size_t m_history;

  static const size_t s_mask = (1 << History) - 1;
};


} // namespace bwtc

#endif
