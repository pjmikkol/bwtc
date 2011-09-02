/**
 * @file DMC.hpp
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
 * Template for Dynamic Markov Chain (DMC) model.
 */


#ifndef BWTC_DMC_HPP
#define BWTC_DMC_HPP

#include "../globaldefs.hpp"
#include "ProbabilityModel.hpp"

#include <vector>
#include <boost/static_assert.hpp>

namespace bwtc {

template <size_t SplitLimit, size_t  MaxStates, size_t Delay, size_t BufferSize,
          typename BitVector, typename BitPredictor>
class DMC : public ProbabilityModel {

 private:
  struct DMCState {
    DMCState()
        : m_zeros(0), m_ones(0), m_zeroTransition(0), m_oneTransition(0) {}
    DMCState(const DMCState& state)
        : m_zeros(state.m_zeros), m_ones(state.m_ones),
          m_zeroTransition(state.m_zeroTransition),
          m_oneTransition(state.m_oneTransition)
    {}

    void delay(DMCState& state) {
      m_zeros /= Delay;
      m_ones /= Delay;
      state.m_zeros = m_zeros;
      state.m_ones = m_ones;
    }
    
    size_t m_zeros;
    size_t m_ones;
    size_t m_zeroTransition;
    size_t m_oneTransition;
    BitPredictor m_predictor;
  };

 public:
  DMC() : m_currentState(0) {
    m_states.reserve(MaxStates);
    m_states.push_back(DMCState());
  }

  ~DMC() {}

  void update(bool bit) {
    if(m_buffer.size() == BufferSize) m_buffer.pop_front();
    m_buffer.push_back(bit);

    DMCState& state = m_states[m_currentState];
    state.m_predictor.update(bit);
    size_t curr;

    // Next state
    if(bit) {
      state = m_states[curr = state.m_oneTransition];
      ++state.m_ones;
    } else {
      state = m_states[curr = state.m_zeroTransition];
      ++state.m_zeros;
    }
    
    if(state.m_ones >= SplitLimit && state.m_zeros >= SplitLimit) {
      size_t n = m_states.size();
      if(n == MaxStates) {
        BitVector oldBuffer(m_buffer);
        resetModel();
        for(size_t i = 0; i < oldBuffer.size(); ++i) {
          update(oldBuffer[i]);
        }
        return;
      } else {
        curr = n;

        m_states.push_back(DMCState(state));
        DMCState& clone = m_states.back();
        clone.delay(state);
        
        if(bit) m_states[m_currentState].m_oneTransition = curr;
        else m_states[m_currentState].m_zeroTransition = curr;
      }
    }
    m_currentState = curr;
  }

  Probability probabilityOfOne() const {
    return m_states[m_currentState].m_predictor.probabilityOfOne();
  }
  
  void resetModel() {
    m_states.clear();
    m_states.reserve(MaxStates);
    m_states.push_back(DMCState());
    m_buffer.clear();
    m_currentState = 0;
  }

  void updateState(bool bit) {
    DMCState& state = m_states[m_currentState];
    if(bit) m_currentState = state.m_oneTransition;
    else m_currentState = state.m_zeroTransition;
  }
  
 private:
  size_t m_currentState;
  BitVector m_buffer;
  std::vector<DMCState> m_states;

  BOOST_STATIC_ASSERT(Delay > 1);
  BOOST_STATIC_ASSERT(BufferSize < SplitLimit*MaxStates - (SplitLimit/Delay)*MaxStates);
};

} //namespace bwtc


#endif
