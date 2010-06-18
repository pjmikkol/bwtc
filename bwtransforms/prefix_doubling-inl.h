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

#ifndef DCSBWT_PREFIX_DOUBLING_INL_H__
#define DCSBWT_PREFIX_DOUBLING_INL_H__

#include <functional>  // for binary_function, bind1st
#include <algorithm>   // for copy, sort, find_if
#include <numeric>     // for partial_sum

#include "prefix_doubling.h"
#include "ternary_partition.h"

namespace dcsbwt {

template <typename StringIterator,
          typename RankIterator,
          typename SuffixIterator>
void SortSuffixesByPrefixDoubling(
    StringIterator str, StringIterator str_end,
    RankIterator rank_array, RankIterator rank_array_end,
    SuffixIterator suffix_array, SuffixIterator suffix_array_end)
{
  uint32 const kMaxLength = 0XfffffffdU;  // == (2^32) - 3
  uint32 length = str_end - str;
  assert(length <= kMaxLength);
  assert(rank_array_end - rank_array == length + 1);
  assert(suffix_array_end - suffix_array == length + 1);
  PrefixDoubler<uint32> pd(length);
  // Construct the initial 1-order
  pd.ComputeOneOrderForCharString(str, str_end, rank_array, rank_array_end,
                                  suffix_array, suffix_array_end);
  // Main repeated doubling procedure
  pd.SortSuffixes(rank_array, rank_array_end,
                  suffix_array, suffix_array_end, 1U);
  // Compute the suffix array from the final rank array
  pd.ComputeFinalSuffixArray(rank_array, rank_array_end,
                             suffix_array, suffix_array_end);
}

namespace {

// Order comparison functor that treats its operands as indices to an array
// and compares the array entries instead of the indices themselves.
template <typename Iterator>
class IndexLess : public std::binary_function<
  typename std::iterator_traits<Iterator>::difference_type,
  typename std::iterator_traits<Iterator>::difference_type, bool>
{
 public:
  typedef typename std::iterator_traits<Iterator>::difference_type Index;
  explicit IndexLess(Iterator begin, Iterator end) : array_(begin) {}
  // The end of the array is not currently stored but a later version
  // could store it and use it for boundary checking (in debug mode).
  bool operator() (Index i, Index j) const { return array_[i] < array_[j]; }
 private:
  Iterator array_;
};

}  // unnamed namespace

template <typename Integer>
template <typename StringIterator, typename RankIterator,
          typename SuffixIterator>
void PrefixDoubler<Integer>::ComputeOneOrder(
    StringIterator str, StringIterator str_end,
    RankIterator rank_array, RankIterator rank_array_end,
    SuffixIterator suffix_array, SuffixIterator suffix_array_end) const
{
  assert(str_end - str == string_length_);
  assert(rank_array_end - rank_array == string_length_ + 1);
  assert(suffix_array_end - suffix_array == string_length_ + 1);
  // Fill with the sequence: string_length_, 0, 1, ..., string_length_-1
  suffix_array[0] = string_length_;  // the empty suffix
  for (SuffixIterator it = suffix_array + 1; it != suffix_array_end; ++it) {
    *it = it - suffix_array - 1;
  }
  // Sort the suffixes into 1-order
  IndexLess<StringIterator> less_for_one_order(str, str_end);
  std::sort(suffix_array + 1, suffix_array_end, less_for_one_order);
  // Compute the rank array
  rank_array[string_length_] = 0;
  SuffixIterator begin, end;
  for (begin = suffix_array + 1; begin != suffix_array_end; begin = end) {
    end = std::find_if(begin, suffix_array_end,
                       std::bind1st(less_for_one_order, *begin));
    // [begin,end) is now an equivalence class
    Integer rank = end - suffix_array - 1;
    for (SuffixIterator it = begin; it != end; ++it) {
      rank_array[*it] = rank;
    }
    if (end - begin == 1) {
      *begin = FinishedSuffixMarker();
    }
    // Guard against infinite loop
    assert(begin - suffix_array < end - suffix_array);
  }
}

template <typename Integer>
template <typename StringIterator, typename RankIterator,
          typename SuffixIterator>
void PrefixDoubler<Integer>::ComputeOneOrderForCharString(
    StringIterator str, StringIterator str_end,
    RankIterator rank_array, RankIterator rank_array_end,
    SuffixIterator suffix_array, SuffixIterator suffix_array_end) const
{
  assert(str_end - str == string_length_);
  assert(rank_array_end - rank_array == string_length_ + 1);
  assert(suffix_array_end - suffix_array == string_length_ + 1);
#ifndef NDEBUG
  // Check that no information is lost by conversion to unsigned char
  typedef typename std::iterator_traits<StringIterator>::value_type CharType;
  for (StringIterator it = str; it != str_end; ++it) {
    assert(*it == static_cast<CharType>(static_cast<unsigned char>(*it)));
  }
#endif
  // Count the number of occurrences of each character
  std::vector<Integer> char_count(256, 0);
  for (StringIterator it = str; it != str_end; ++it) {
    unsigned int ch = *it;
    ++char_count[ch];
  }
  // Compute the rank for suffixes starting with each character
  std::vector<Integer> rank(char_count);
  std::partial_sum(char_count.begin(), char_count.end(),
                   rank.begin());
  assert(rank[255] == string_length_);
  // The rank is also the suffix array position of
  // the last suffix in the equivalence class.
  // A shift by one gives us the last position of the previous
  // equivalence class.
  std::vector<uint32> suffix_array_position(256);
  suffix_array_position[0] = 0;
  std::copy(rank.begin(), rank.end() - 1,
            suffix_array_position.begin() + 1);
  for (StringIterator it = str; it != str_end; ++it) {
    unsigned int ch = *it;
    Integer str_pos = it - str;
    rank_array[str_pos] = rank[ch];
    suffix_array[++suffix_array_position[ch]] = str_pos;
  }
  assert(rank == suffix_array_position);  // !!
  // the empty suffix
  rank_array[string_length_] = 0;
  suffix_array[0] = string_length_;
  // Could add finished suffix markers here but is it worth the effort?
}

template <typename Integer>
template <typename RankIterator, typename SuffixIterator>
void PrefixDoubler<Integer>::ComputeInitialSuffixArray(
    RankIterator rank_array, RankIterator rank_array_end,
    SuffixIterator suffix_array, SuffixIterator suffix_array_end,
    Integer order_depth) const
{
  // Fill with the sequence: string_length_, 0, 1, ..., string_length_-1
  suffix_array[0] = string_length_;  // the empty suffix
  for (SuffixIterator it = suffix_array + 1; it != suffix_array_end; ++it) {
    *it = it - suffix_array - 1;
  }
  RefineSubrange(rank_array, suffix_array,
              suffix_array + 1, suffix_array_end, 0);
}

template <typename Integer>
template <typename RankIterator, typename SuffixIterator>
void PrefixDoubler<Integer>::ComputeFinalSuffixArray(
    RankIterator rank_array, RankIterator rank_array_end,
    SuffixIterator suffix_array, SuffixIterator suffix_array_end) const
{
#ifndef NDEBUG
  // markers identify uninitialized entries
  std::fill(suffix_array, suffix_array_end, FinishedSuffixMarker());
#endif
  for (RankIterator it = rank_array; it != rank_array_end; ++it) {
    suffix_array[*it] = it - rank_array;
  }
#ifndef NDEBUG
  SuffixIterator uninitialized_entry =
      std::find(suffix_array, suffix_array_end, FinishedSuffixMarker());
  assert(suffix_array_end - uninitialized_entry == 0);
#endif
}

// Refine the order of a group of consecutive suffixes.
// This is where the actual modification of ranks takes place.
// Preconditions (situation on entry):
// 1. [begin,end) is a sub range of the suffix array.
//    The suffixes in this subrange have a common prefix of length
//    common_prefix_length. The common_prefix_length can be 0.
// 2. The subrange occupies its final place in the suffix array,
//    that is, the final suffix array positions of the suffixes
//    are in this subrange. More formally, the number of other suffixes
//    that are smaller than all the suffixes in this subrange is
//    begin-suffix_array, and the number of other suffixes
//    that are larger than all the suffixes in this subrange is
//    suffix_array_end-end.
// 3. The other suffixes that are smaller than the suffixes in this
//    subrange have ranks smaller than begin-suffix_array.
//    The other suffixes that are larger than the suffixes in this
//    subrange have ranks larger than or equal to suffix_array_end-end.
//    (Together the above conditions imply that the ranks of the suffixes
//    in this subrange are larger than or equal to begin-suffix_array
//    and smaller than or equal to the smallest rank of a larger suffix.)
// 4. The rank array represents some weak order consistent with
//    the lexicographic order, and thus a refinement of some h-order.
//    The value of h is not provided (or needed), but it is assumed
//    that the caller knows such h.
//    Often h==common_prefix_length but this is not required.
// All the preconditions hold on exit, too, as well as the following
// additional condition:
// * The ranks of the suffixes in the subrange are now consistent
//   with the (common_prefix_length+h)-order and their order in the
//   subrange is consistent with the ranks.
template <typename Integer>
template <typename RankIterator, typename SuffixIterator>
void PrefixDoubler<Integer>::RefineSubrange(
    RankIterator rank_array, SuffixIterator suffix_array,
    SuffixIterator begin, SuffixIterator end,
    Integer common_prefix_length) const
{
  static const int kSmallRangeLimit = 15;
  assert(end - begin >= 0);
  RankAtDepth<RankIterator> sortkey(rank_array, common_prefix_length);
  while (end - begin > kSmallRangeLimit) {
    // Divide the range into three subranges:
    // smaller-than-pivot, equal-to-pivot, and larger-than-pivot
    SuffixIterator pivot = ChoosePivot(begin, end, sortkey);
    std::pair<SuffixIterator, SuffixIterator> equal_group =
        TernaryPartition(begin, end, pivot, sortkey);
    SuffixIterator equal_begin = equal_group.first;
    SuffixIterator equal_end = equal_group.second;
    // Partition the smaller-than-pivot range recursively.
    RefineSubrange(rank_array, suffix_array, begin, equal_begin,
                   common_prefix_length);
    // Assign a new rank for the equal-to-pivot group
    Integer new_rank = equal_end - suffix_array - 1;
    for (SuffixIterator it = equal_begin; it != equal_end; ++it) {
      rank_array[*it] = new_rank;
    }
    if (equal_end - equal_begin == 1) {
      *(equal_begin) = FinishedSuffixMarker();
    }
    // Loop (tail recurse) to partition the larger-than-pivot group.
    begin = equal_end;
    assert(end - begin >= 0);
  }
  if (begin == end) return;
  if (end - begin == 1) {
    rank_array[*begin] = begin - suffix_array;
    *begin = FinishedSuffixMarker();
    return;
  }
  RefineSmallSubrange(rank_array, suffix_array,
                      begin, end, common_prefix_length);
}

template <typename Integer>
template <typename RankIterator, typename SuffixIterator>
void PrefixDoubler<Integer>::RefineSmallSubrange(
    RankIterator rank_array, SuffixIterator suffix_array,
    SuffixIterator begin, SuffixIterator end,
    Integer common_prefix_length) const
{
  // Selection sort for a small range
  assert(end - begin >= 0);
  RankAtDepth<RankIterator> sortkey(rank_array, common_prefix_length);
  SuffixIterator smallest = begin;
  Integer smallest_key = sortkey(*smallest);
  for (SuffixIterator it = begin + 1; it < end; ++it) {
    if (sortkey(*it) < smallest_key) {
      smallest_key = sortkey(*it);
      smallest = it;
    }
  }
  using std::swap;
  swap(*smallest, *begin);
  SuffixIterator previous = begin;
  Integer previous_key = smallest_key;
  for (SuffixIterator it = begin + 1; it < end; ++it) {
    smallest = it;
    smallest_key = sortkey(*it);
    for (SuffixIterator it2 = it + 1; it2 < end; ++it2) {
      if (sortkey(*it2) < smallest_key) {
        smallest_key = sortkey(*it2);
        smallest = it2;
      }
    }
    swap(*smallest, *it);
    if (smallest_key > previous_key) {
      Integer rank = (it - 1) - suffix_array;
      for (SuffixIterator it2 = previous; it2 < it; ++it2) {
        rank_array[*it2] = rank;
      }
      if (it - previous == 1) {
        *previous = FinishedSuffixMarker();
      }
      previous = it;
      previous_key = smallest_key;
    }
  }
  Integer rank = (end - 1) - suffix_array;
  for (SuffixIterator it = previous; it < end; ++it) {
    rank_array[*it] = rank;
  }
  if (end - previous == 1) {
    *previous = FinishedSuffixMarker();
  }
}

// Go from h-order to 2h-order for all suffixes.
// This is done by finding all equivalence classes that have not been
// marked finished and calling DoubleGroup() for them.
// Consecutive marked suffixes or groups are combined.
// (DoubleGroup() does the initial marking of finished suffixes.)
template <typename Integer>
template <typename RankIterator, typename SuffixIterator>
Integer PrefixDoubler<Integer>::Double(
    RankIterator rank_array, RankIterator rank_array_end,
    SuffixIterator suffix_array,  SuffixIterator suffix_array_end,
    Integer order_depth) const
{
  assert(order_depth > 0);
  Integer size = suffix_array_end - suffix_array;
  assert(size == string_length_ + 1);
  assert(rank_array_end - rank_array == size);
  bool previous_is_finished = false;
  Integer previous_finished_first = 0;
  Integer first = 0;  // first position in equivalence class or finished group
  Integer next = size;  // first position in the next class or group
  for (first = 0; first < size; first = next) {
    assert(first <= string_length_);
    assert(first >= 0);
    Integer suffix = suffix_array[first];
    if (IsFinishedSuffixOrGroup(suffix)) {
      if (IsFinishedGroup(suffix)) {
        next = suffix_array[first + 1];
      } else {
        assert(IsFinishedSuffix(suffix));
        next = first + 1;
      }
      if (previous_is_finished) {
        // Combine with a preceding finished suffix or group
        suffix_array[previous_finished_first] = FinishedGroupMarker();
        suffix_array[previous_finished_first + 1] = next;
      } else {
        previous_is_finished = true;
        previous_finished_first = first;
      }
    } else {
      assert(suffix <= string_length_);
      assert(suffix >= 0);
      Integer rank = rank_array[suffix];
      assert(rank <= string_length_);
      assert(rank >= 0);
      next = rank + 1;  // rank is the position of the last suffix
                        // in the equivalence class
      RefineSubrange(rank_array, suffix_array,
                  suffix_array + first, suffix_array + next, order_depth);
      previous_is_finished = false;
    }
    // guard against infinite loop
    assert(next > first);
  }
  return order_depth * 2;
}

// Sort all suffixes by repeated doubling.
template <typename Integer>
template <typename RankIterator, typename SuffixIterator>
void PrefixDoubler<Integer>::SortSuffixes(
    RankIterator rank_array, RankIterator rank_array_end,
    SuffixIterator suffix_array, SuffixIterator suffix_array_end,
    Integer order_depth) const
{
  assert(order_depth > 0);
  while (order_depth < string_length_) {
    order_depth = Double(rank_array, rank_array_end,
                         suffix_array, suffix_array_end, order_depth);
  }
}

} // namespace dcsbwt

#endif  // DCSBWT_PREFIX_DOUBLING_INL_H__
