/**
 * @file longsequences.cc
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
 * Header for detecting and replacing long repetitive sequences.
 */

#ifndef BWTC_LONGSEQUENCES_H_
#define BWTC_LONGSEQUENCES_H_

#include <vector>
#include "../globaldefs.hpp"

namespace bwtc {

namespace long_sequences {

/**
 * chunk is used to represent single (relatively small) fixed size block of
 * input sequence.
 *
 * This struct is used by SequenceDetector- and Expander-objects.
 */
struct chunk {
  chunk() : position(0), hash_value(0) {}
  chunk(uint32 p, uint32 h) : position(p), hash_value(h) {}
  uint32 position;
  uint32 hash_value;
};

struct cmp_chunk {
  cmp_chunk(uint32 *c) : counts(c) {}
  bool operator()(const chunk& a, const chunk& b) {
    if (counts[a.hash_value] < counts[b.hash_value])
      return true;
    else if (counts[a.hash_value] > counts[b.hash_value])
      return false;
    else if (a.hash_value < b.hash_value)
      return true;
    else if (a.hash_value > b.hash_value)
      return false;
    else return a.position < b.position;
  }
  uint32 *counts;
};


/** Calculates frequencies from range [source_begin, source_end)
 *
 *  @param target where the frequency information is stored
 *  @param source_begin start of the inspected sequence
 *  @param source_begin one element past the end of the inspected sequence
 */
template <typename T>
void CalculateFrequencies(T *target, byte *source_begin, byte *source_end) {
  while(source_begin != source_end) ++target[*source_begin++];
}

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
