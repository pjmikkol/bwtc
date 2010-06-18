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

// Prefix doubling is a technique for sorting the suffixes
// of a string of length n in O(n log n) worst-case time.
// This file contains an implementation of a prefix doubling algorithm
// based on the Larsson-Sadakane algorithm described in
//
//   J. Larsson, K. Sadakane: Faster Suffix Sorting.
//   Theoretical Computer Science, in press 2007.
//   http://eprints.csce.kyushu-u.ac.jp/33/
//
// The simple way to use the algorithm is the function
// SortSuffixesByPrefixDoubling().
// More flexibility and more chances for optimizations are offered
// by the class PrefixDoubler.

#ifndef DCSBWT_PREFIX_DOUBLING_H__
#define DCSBWT_PREFIX_DOUBLING_H__

#include <functional>  // for unary_function
#include <cassert>

namespace dcsbwt {

// The function SortSuffixesByPrefixDoubling takes three arrays as
// an argument. The first array contains the input string, and the
// other two contain the output on exit:
// - The rank array contains a rank (an integer in the range [0..n])
//   for each suffix such that the order of the ranks is the same as
//   the lexicographic order of the suffixes.
// - The suffix array contains the suffixes in the lexicographic order
//   with each suffix represented by its starting position.
// For example:
//   input string: BANANA
//     rank array: 4362510
//   suffix array: 6531042
//
// The arrays are, in fact, STL random access iterator ranges.
// Here is a precise description of the input and the output:
// - The range [str,str_end) contains the input string; it will not be
//   modified. The length of the range, i.e., the length of the string,
//   (call it n) must be at most (2^32)-3.
//   The type of the elements should be either char or unsigned char.
//   (Other types are possible as long as they can be converted to
//   unsigned char; this is checked in debug mode.) The characters
//   are converted to unsigned char for purposes of comparison.
//   (The standard string comparison operators and strcmp-like functions
//   behave the same way.)
// - The range [rank_array,rank_array_end) contains the ranks on output;
//   the contents on input are not used. The size of the range should
//   be n+1 (one larger than the size of [str,str_end)). The element
//   type should be an integral type capable of storing all values in
//   the range [0..n] and supporting conversion to and from uint32.
//   After execution, rank_array[i] contains the number of the suffixes
//   of the string (including the empty suffix and the full string)
//   that are smaller than the suffix [str+i,str_end).
//   An alternative definition is that the range contains a permutation
//   of the integers [0..n] satisfying for all i,j in [0..n]:
//     rank_array[i] < rank_array[j] if and only if
//     std::lexicographical_compare(str+i, str_end, str+j, str_end,
//                                  std::less<unsigned char>())
//     returns true.
// - The range [suffix_array,suffix_array_end) contains the suffix array
//   on output; the contents on input are not used. The size of the range
//   should be n+1 (the same as [rank_array,rank_array_end)). The element
//   type should be an integral type capable of storing all values in
//   the range [0..n+2] (the values n+1 and n+2 are needed by internal
//   computation) and supporting conversion to and from uint32.
//   After execution, the range contains a permutation of
//   the integers [0..n], where a value i represents the suffix
//   [str+i,str_end) and the order is the same as the lexicographical
//   order of the suffixes. It satisfies, for all i in [0..n]:
//   rank_array[suffix_array[i]] == suffix_array[rank_array[i]] == i.
//   (If you wonder why both rank and suffix array are produced when
//   one can so easily be constructed from the other, it is because both
//   arrays are needed during the computation, so one can as well
//   use them both for output.)

template <typename StringIterator,
          typename RankIterator,
          typename SuffixIterator>
void SortSuffixesByPrefixDoubling(
    StringIterator str, StringIterator str_end,
    RankIterator rank_array, RankIterator rank_array_end,
    SuffixIterator suffix_array, SuffixIterator suffix_array_end);

// SortSuffixesByPrefixDoubling is not sufficiently flexible for some
// special situations (in particular, those arising in the suffix sorting
// algorithms based on difference covers). The class PrefixDoubler provides
// the core algorithm in its bare form together with some supporting
// components that can be used for building special applications.
// To understand PrefixDoubler you need to know more about the algorithm.
//
// An h-order of the suffixes is an order based on the first h characters
// (call it the h-prefix) of each suffix. If the suffix is shorter than h,
// the h-prefix is the suffix itself.  An h-order is not a total order
// in general (but a weak order, see below) as multiple suffixes can have
// the same h-prefix. For example, under 0-order all suffixes are equal.
// However, for h >= n, the h-order is always the lexicographical total
// order of the suffixes.
//
// The algorithm is based on a doubling procedure that turns an h-order
// into a 2h-order. The outline of the algorithm is:
// 1. Obtain an initial h-order for some h > 0.
// 2. Apply the doubling procedure repeatedly until fully sorted.
// The second step is the core algorithm implemented by the member function
// PrefixDoubler::SortSuffixes(). It takes as an input the value of h (called
// order_depth in the code) and the two arrays, rank and suffix array,
// which must be consistent with the initial h-order on input.
// On output, the rank array is as with SortSuffixesByPrefixDoubling()
// but the suffix array is not (see below). The suffix array can be
// computed from the rank array using PrefixDoubler::ComputeFinalSuffixArray().
// Note that the original string is not an input to the main algorithm;
// it is only used for obtaining the initial h-order. There are functions
// that help in computing the initial h-order, but users can also write
// their own, which is the main source for the added flexibility.
//
// The rank array is the primary representation of the h-order at all stages
// of the algorithm. It must satisfy the following conditions:
// 1. If rank_array[i] < rank_array[j], then the suffix starting at
//    position i is lexicographically smaller than the suffix starting at
//    position j.
// 2. If rank_array[i] == rank_array[j], then suffixes i and j have
//    the same h-prefix. (The reverse implication is not required;
//    see the paragraph on weak orders below.)
// 3. rank_array[i] is equal to the number of other suffixes whose rank is
//    smaller or equal to rank_array[i]. For example, rank_array[n]=0
//    for all h > 0 (n is the empty suffix).
// When all suffixes have a different rank (e.g. for all h >= n), this
// definition coincides with the earlier definition of the output
// of SortSuffixesByPrefixDoubling().
//
// Example:  input string: BANANA
//     1-order rank array: 4363630  <-- initial 1-order
//     2-order rank array: 4363610  <-- after first doubling
//     4-order rank array: 4362510  <-- after second doubling
//                                      this is the total order
//
// The suffix array in its basic form contains the suffixes (i.e.,
// a permutation of [0..n]) in the order defined by the rank array,
// with equal suffixes in arbitrary order.
//
// Example:  input string: BANANA
//     1-order rank array: 4363630
//   1-order suffix array: 6135024  <-- both of these would be valid
//                         6531042  <-- inputs to SortSuffixes()
//
// An advanced form allows overwriting some entries with special markers
// as explained later. For basic use it is enough to know that:
// - On input to the main algorithm SortSuffixes(), the suffix array
//   can but does not have to contain the markers.
// - On output, the suffix array contains markers to a large extent,
//   and the funtion ComputeFinalSuffixArray() should be used for
//   obtaining the final suffix array if needed.
//
// The class PrefixDoubler stores no information about the state of the
// algorithm except the length of the string; each function takes all
// necessary information as an argument. This design gives the user full
// control over the storage and representation of the arrays, which can
// be useful since these arrays can be really large. The reasons
// why this class exists at all (instead of separate functions) are:
// - collecting the related functions together
// - encapsulation of the special markers explained below
// - parametrization of the integral type
// The last point refers to the template paraneter type Integer, which is
// used in all internal integral computations. It should be able to represent
// all values in the range [0..n+2]. The value type of the rank and
// suffix arrays should support conversion to and from the type Integer
// for all values in the range [0..n+2].
//
// The rest of the documentation is not essential for basic use but
// may be useful for understanding the code and for some optimizations.
//
// A weak order over a set is an equivalence class partition of the set
// with a total order on the equivalence classes. Every h-order is
// a weak order over the set of suffixes. Weak orders are exactly
// those that fit the STL requirement of strict weak ordering (see 25.3
// in the C++ Standard). Two relations on weak orders are of interest to us:
// - Two weak orders W1 and W2 are *consistent* if there are no two
//   elements x and y for which x<y according to W1 and y<x according
//   to W2. All h-orders are consistent with each other.
// - A weak order W1 is a *refinement* of a weak order W2 if W1 and W2
//   are consistent and the equivalence class partition of W1 is a
//   refinement of the equivalence class partition of W2.
//   A g-order is a refinement of an h-order if g > h.
// The rank array is able to represent any weak order on the suffixes
// A valid order for a specific value of h is any order that:
// - is consistent with the lexicographical total order and
// - is a refinement of the h-order.
// These are the conditions 1 and 2 in the earlier rank array definition.
// Any valid order for h > 0 is a legal input to the algorithm.
// A more refined order makes the algorithm faster in general.
// Note though that an additional refinement is much more effective if
// it allows increasing the value of h. Additional refinements have to
// be reflected in the suffix array, too.
//
// Example:  input string: BANANA
//     1-order rank array: 4363630  <-- no additional refinement
//                         4363610  <-- additional refinement
//   1-order suffix array: 6531042  <-- valid for both cases
//                         6135024  <-- not valid for refined case
//
// A suffix that has been separated from other suffixes (i.e.,has
// a unique rank), is called finished: its rank and its position in the
// suffix array do not change anymore. The suffix array entry for
// a finished suffix is not needed anymore (the rank is enough)
// and may be overwritten by special a special marker,
// whose purpose is to quickly jump over the finished parts.
// Two kinds of markings are used:
// - The value FinishedSuffixMarker() marks a single finished suffix.
// - The value FinishedGroupMarker() marks the beginning of a group of
//   more than one consecutive finished suffixes. The entry following
//   the marker contains the next position after the end of the group.
//   The remaining entries in the group are never accessed.
// Any combination of unmarked finished suffixes, single marked
// suffixes and marked groups are legal including consecutive marked
// suffixes and groups. Unmarked finished suffixes get marked and
// consecutive marked suffixes and groups get combined when encountered.
// Adding markers already during the computation of the initial
// h-order may speed up the algorithm. (Besides allowing jumping over
// groups, the markers avoid cache misses due to rank array accesses.)

template <typename Integer>
class PrefixDoubler {
 public:
  // Only the length of the string is set on initialization
  explicit PrefixDoubler(Integer string_length)
      : string_length_(string_length) {}
  // no default constructor
  // implicit destructor
  // implicit copy constructor
  // implicit copy assignment

  // The main repeated doubling algorithm
  template <typename RankIterator, typename SuffixIterator>
  void SortSuffixes(
      RankIterator rank_array, RankIterator rank_array_end,
      SuffixIterator suffix_array, SuffixIterator suffix_array_end,
      Integer order_depth) const;

  // Compute an initial 1-order from a string.
  // This works for any string whose characters can be compared.
  // The ordering is based on std::less<value_type> for the
  // value_type of the StringIterator.
  // WARNING: If the value type is char and the string contains values
  // outside [0..128), the order may be inconsistent with the standard
  // string comparisons and strcmp-like functions, which convert to
  // unsigned char for comparisons. Use ComputeOneOrderForCharString()
  // instead.
  // TODO: add support for user-defined character comparison operator
  template <typename StringIterator, typename RankIterator,
            typename SuffixIterator>
  void ComputeOneOrder(
      StringIterator str, StringIterator str_end,
      RankIterator rank_array, RankIterator rank_array_end,
      SuffixIterator suffix_array, SuffixIterator suffix_array_end) const;

  // Compute an initial 1-order from a string of (unsigned) char.
  // This works for strings whose characters are convertible to unsigned char.
  // In debug mode, there is a check that the conversion loses no information.
  // The ordering is based on std::less<unsigned char>.
  // This is probably faster than ComputeOneOrder(), and produces an order
  // consistent with the standard string comparisons and strcmp-like functions.
  template <typename StringIterator, typename RankIterator,
            typename SuffixIterator>
  void ComputeOneOrderForCharString(
      StringIterator str, StringIterator str_end,
      RankIterator rank_array, RankIterator rank_array_end,
      SuffixIterator suffix_array, SuffixIterator suffix_array_end) const;

  // Compute the initial suffix_array including markers
  // from the initial rank array.
  template <typename RankIterator, typename SuffixIterator>
  void ComputeInitialSuffixArray(
      RankIterator rank_array, RankIterator rank_array_end,
      SuffixIterator suffix_array, SuffixIterator suffix_array_end,
      Integer order_depth) const;

  // Compute the final suffix array (without markers) from the final
  // rank array.
  // This should only be called with the fully sorted rank array
  // (checked in debug mode).
  template <typename RankIterator, typename SuffixIterator>
  void ComputeFinalSuffixArray(
      RankIterator rank_array, RankIterator rank_array_end,
      SuffixIterator suffix_array, SuffixIterator suffix_array_end) const;

  // Special markers for finished suffixes
  inline Integer FinishedSuffixMarker() const { return string_length_ + 1; }
  inline Integer FinishedGroupMarker() const { return string_length_ + 2; }
  inline bool IsFinishedSuffix(Integer i) const {
    return FinishedSuffixMarker() == i;
  }
  inline bool IsFinishedGroup(Integer i) const  {
    return FinishedGroupMarker() == i;
  }
  inline bool IsFinishedSuffixOrGroup(Integer i) const  {
    assert(i <= string_length_ + 2);
    return i > string_length_;
  }

  // Do one doubling step.
  // Returns the new order_depth.
  template <typename RankIterator, typename SuffixIterator>
  Integer Double(
      RankIterator rank_array, RankIterator rank_array_end,
      SuffixIterator suffix_array, SuffixIterator suffix_array_end,
      Integer order_depth) const;

 private:
  Integer string_length_;

  // Functor that retrieves the rank at i+order_depth
  template <typename RankIterator>
  class RankAtDepth : public std::unary_function<Integer, Integer> {
   public:
    RankAtDepth(RankIterator rank_array, Integer order_depth)
        : rank_plus_depth_(rank_array + order_depth) {}
    Integer operator() (Integer i) const {
      return rank_plus_depth_[i];
    }
   private:
    RankIterator rank_plus_depth_;
  };

  // Do one doubling step for a group of suffixes.
  template <typename RankIterator, typename SuffixIterator>
  void RefineSubrange(
      RankIterator rank_array, SuffixIterator suffix_array,
      SuffixIterator begin, SuffixIterator end,
      Integer common_prefix_length) const;

  // Do one doubling step for a group of suffixes.
  template <typename RankIterator, typename SuffixIterator>
  void RefineSmallSubrange(
      RankIterator rank_array, SuffixIterator suffix_array,
      SuffixIterator begin, SuffixIterator end,
      Integer common_prefix_length) const;
};

}  // namespace dcsbwt

#endif  // DCSBWT_PREFIX_DOUBLING_H__
