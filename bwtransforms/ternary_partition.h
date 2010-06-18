// Copyright 2007 Google Inc.

// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

// Ternary partition algorithm. Includes main functions:
//   ChoosePivot()
//   TernaryPartition()
// and supporting code
//   class RandomNumberGenerator
//   function MedianOfThree()
//
// Genereric comparison is achieved using a functor that returns
// the comparison key instead of a comparison functor.

#ifndef DCSBWT_TERNARY_PARTITION_H__
#define DCSBWT_TERNARY_PARTITION_H__

#include "inttypes.h"

#include <utility>    // for pair, make_pair
#include <algorithm>  // for swap
#include <cstdlib>    // for drand48

namespace dcsbwt {

namespace ternary_partition {

// Currently implemented using drand48.
// TODO: Switch to an implementation using Mersenne Twister or
// to using directly TR1 random number generators.
class RandomNumberGenerator {
 public:
  RandomNumberGenerator() { srand48_r(87989823L, &buffer_); }
  typedef uint32 result_type;
  static const result_type min = 0;
  static const result_type max = (1UL << 31) - 1;
  result_type operator()() {
    long int result = 0;
    lrand48_r(&buffer_, &result);
    return result;
  }
  result_type Uniform(uint32 range) {
    assert(range - 1 <= max);
    return operator()() % range;
  }
 private:
  drand48_data buffer_;
};

template <typename Iterator, typename Key>
Iterator MedianOfThree(Iterator a, Iterator b, Iterator c, Key key) {
  typename Key::result_type akey = key(*a);
  typename Key::result_type bkey = key(*b);
  typename Key::result_type ckey = key(*c);
  if (akey < bkey) {
    if (bkey < ckey) return b;
    else if (akey < ckey) return c;
    else return a;
  } else {
    if (akey < ckey) return a;
    else if (bkey < ckey) return c;
    else return b;
  }
}

}  // namespace ternary_partition

// Choose a pivot for partitioning.
// Returns an iterator to the chosen pivot.
template <typename Iterator, typename Key>
Iterator ChoosePivot(Iterator begin, Iterator end, Key key) {
  using ternary_partition::MedianOfThree;
  static ternary_partition::RandomNumberGenerator rng;
  int32 size = end - begin;
  if (size < 100) {
    return begin + rng.Uniform(size);
  //return begin + rng.Uniform(std::distance(begin,end));
  } else if (size < 1000) {
    return MedianOfThree(begin + rng.Uniform(size),
                         begin + rng.Uniform(size),
                         begin + rng.Uniform(size),
                         key);
  } else {
    // pseudo-median of nine
    return MedianOfThree(
        MedianOfThree(begin + rng.Uniform(size),
                      begin + rng.Uniform(size),
                      begin + rng.Uniform(size),
                      key),
        MedianOfThree(begin + rng.Uniform(size),
                      begin + rng.Uniform(size),
                      begin + rng.Uniform(size),
                      key),
        MedianOfThree(begin + rng.Uniform(size),
                      begin + rng.Uniform(size),
                      begin + rng.Uniform(size),
                      key),
        key);
  }
}

// Permute the elements in [begin,end) so that elements smaller than
// pivot are on the left, elements equal to pivot are in the middle,
// and elements larger than pivot are on the right.
// Returns a range identifying the equal part.
template <typename Iterator, typename Key>
std::pair<Iterator,Iterator> TernaryPartition(Iterator begin, Iterator end,
                                              Iterator pivot, Key key) {
  // The unqualified call of swap (i.e. without std::) ensures
  // that overloads of swap for user-defined types are found.
  // This using declaration ensures that the spezializations
  // for built-in types in std are found as well.
  using std::swap;
  typename Key::result_type pivot_key = key(*pivot);
  Iterator smaller_end = begin;
  Iterator larger_begin = end;
  Iterator current = begin;
  while (current != larger_begin) {
    typename Key::result_type current_key = key(*current);
    if (current_key < pivot_key) {
      swap(*current++, *smaller_end++);
    } else if (pivot_key < current_key) {
      swap(*current, *--larger_begin);
    } else {
      ++current;
    }
  }
  return std::make_pair(smaller_end, larger_begin);
}

}  // namespace dcsbwt

#endif  // DCSBWT_TERNARY_PARTITION_H__
