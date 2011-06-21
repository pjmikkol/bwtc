/**
 * @file ProbabilityModel.cpp
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
 * Base class for probability models.
 */

#include <algorithm> // for std::fill
#include <iostream>

#include "../globaldefs.hpp"
#include "ProbabilityModel.hpp"
#include "BitPredictors.hpp"
#include "FSM.hpp"

namespace bwtc {

ProbabilityModel* GiveProbabilityModel(char choice) {
  switch(choice) {
    case 'm':
      if( verbosity > 1)
        std::clog << "Remembering 8 previous bits\n";
      return new SimpleMarkov<byte>();
    case 'M':
      if( verbosity > 1)
        std::clog << "Remembering 16 previous bits\n";
      return new SimpleMarkov<unsigned short int>();      
    case 'b':
      if( verbosity > 1)
        std::clog << "Remembering 4 previous bits.\n";
      return new BitPredictor<4>();
    case 'B':
      if( verbosity > 1)
        std::clog << "Using FSM.\n";
      return new FSM<6, BitPredictor<8> >();
    case 'n':
    default:
      if( verbosity > 1)
        std::clog << "Remembering single previous bit\n";
      return new ProbabilityModel();
  }
}


/*************************************************************************
 * SimpleMarkov: An example of how to integrate new probability model    *
 * to program.                                                           *
 *                                                                       * 
 * Simple template-based probability-model which remembers               *
 * 8*sizeof(Integer) previous bits. Consumes huge amount of memory       *
 * 2^(8*sizeof(Integer) bytes) so practically this is  usable  only with *
 * bytes and short integers.                                             *
 *************************************************************************/
template <typename UnsignedInt>
SimpleMarkov<UnsignedInt>::SimpleMarkov() :
    m_prev(static_cast<UnsignedInt>(0)), m_history(0)
{
  uint64 size = (static_cast<uint64>(1) << 8*sizeof(UnsignedInt)) - 1;
  m_history = new char[size];
  std::fill(m_history, m_history + size, 0);
}

template <typename UnsignedInt>
SimpleMarkov<UnsignedInt>::~SimpleMarkov() {
  delete [] m_history;
}

template <typename UnsignedInt>
void SimpleMarkov<UnsignedInt>::update(bool bit) {
  if (bit) {
    if (m_history[m_prev] < 2 )
      ++m_history[m_prev];
  }
  else {
    if(m_history[m_prev] > -2)
      --m_history[m_prev];
  }
  m_prev <<= 1;
  m_prev |= (bit)? 1 : 0;
}

template <typename UnsignedInt>
Probability SimpleMarkov<UnsignedInt>::probabilityOfOne() const {
  Probability val = kProbabilityScale >> (kLogProbabilityScale/2);
  if (m_history[m_prev] > 0) return val << 2*m_history[m_prev];
  else return val >> 2*m_history[m_prev];
}

template <typename UnsignedInt>
void SimpleMarkov<UnsignedInt>::resetModel() {
  /* Seems to work better when not resetting the model for different
   * contexts. */
  uint64 size = (static_cast<uint64>(1) << 8*sizeof(UnsignedInt)) - 1;
  std::fill(m_history, m_history + size, 0);
}

} //namespace bwtc
