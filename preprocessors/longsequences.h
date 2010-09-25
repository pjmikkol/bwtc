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

#ifndef BWTC_LONGSEQUENCES_H_
#define BWTC_LONGSEQUENCES_H_

#include <vector>

#include "../globaldefs.h"

namespace bwtc {

namespace long_sequences {

/* To avoid costly memory-allocations we allocate the array used in *
 * computations of Borders only once.                               */
template <typename T>
class Border {
 public:
  explicit Border(unsigned length) : length_(length) {
    border_ = new int[length_ + 1];
  }

  ~Border() {
    delete [] border_;
  }

  int operator()(const T *source) {
    border_[0] = -1;
    int t = -1;
    for(unsigned j = 1; j <= length_; ++j) {
      while(t >= 0 && source[t] != source[j-1]) t = border_[t];
      border_[j] = ++t;
    }
    return border_[length_];
  }

 private:
  int *border_;
  unsigned length_;
};

} //namespace long_sequences

uint64 CompressSequences(byte *from, uint64 length, int memory_constraint,
                         unsigned window_size, int threshold);

} 

#endif
