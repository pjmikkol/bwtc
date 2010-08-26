/**************************************************************************
 *  Copyright 2010, Pekka Mikkola, pjmikkol (at) cs.helsinki.fi           *
 *                                                                        *
 *  This file is part of bwtc.                                            *
 *                                                                        *
 *  bwtc is free software: you can redistribute it and/or modify          *
 *  it under the terms of the GNU General Public License as published by  *
 *  the Free Software Foundation, either version 3 of the License, or     *
 *  (at your option) any later version.                                   *
 *                                                                        *
 *  bwtc is distributed in the hope that it will be useful,               *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *  GNU General Public License for more details.                          *
 *                                                                        *
 *  You should have received a copy of the GNU General Public License     *
 *  along with bwtc.  If not, see <http://www.gnu.org/licenses/>.         *
 **************************************************************************/

#include <algorithm> // for std::fill
#include <iostream>

#include "../globaldefs.h"
#include "base_prob_model.h"

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
    prev_(static_cast<UnsignedInt>(0)), history_(NULL)
{
  uint64 size = (static_cast<uint64>(1) << 8*sizeof(UnsignedInt)) - 1;
  history_ = new char[size];
  std::fill(history_, history_ + size, 0);
}

template <typename UnsignedInt>
SimpleMarkov<UnsignedInt>::~SimpleMarkov() {
  delete [] history_;
}

template <typename UnsignedInt>
void SimpleMarkov<UnsignedInt>::Update(bool bit) {
  if (bit) {
    if (history_[prev_] < 2 )
      ++history_[prev_];
  }
  else {
    if(history_[prev_] > -2)
      --history_[prev_];
  }
  prev_ <<= 1;
  prev_ |= (bit)? 1 : 0;
}

template <typename UnsignedInt>
Probability SimpleMarkov<UnsignedInt>::ProbabilityOfOne() {
  Probability val = kProbabilityScale >> (kLogProbabilityScale/2);
  if (history_[prev_] > 0) return val << 2*history_[prev_];
  else return val >> 2*history_[prev_];
}

template <typename UnsignedInt>
void SimpleMarkov<UnsignedInt>::ResetModel() {
  /* Seems to work better when not resetting the model for different
   * contexts. */
  uint64 size = (static_cast<uint64>(1) << 8*sizeof(UnsignedInt)) - 1;
  std::fill(history_, history_ + size, 0);
}

} //namespace bwtc
