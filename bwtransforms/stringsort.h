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

#ifndef DCSBWT_STRINGSORT_H__
#define DCSBWT_STRINGSORT_H__

#include "inttypes.h"

namespace dcsbwt {

////////////////////////////////////////////////////////////////
// Sort a set of suffixes of a text using string sorting algorithms,
// i.e., algorithms that do not take advantage of the fact the input
// strings are suffixes of a single string. The primary algorithm
// is multikey quicksort.
//
// ARGUMENTS
//
// [text,text_end) contains the text whose suffixes are being sorted
//   CharIterator must be a random access iterator with a value type
//   convertible to unsigned char.
//
// [suffix_area,suffix_area_end) contains the suffixes and possibly
//   additional working space. A value i represents [text+i,text_end).
//   UnsignedInteger must be a built-in unsigned integer type.
//   * On entry, the suffixes to be sorted must be at the end of
//     the suffix area in [suffixes_in,suffix_area_end).
//     The contents of [suffix_area,suffixes_in) are not used.
//   * On exit, the suffixes are at the beginning of the suffix area,
//     i.e, in [suffix_area,suffixes_end), where
//     suffixes_end = suffix_area + (suffix_area_end - suffixes_in)
//     is the return value of the function.
//     The contents of the rest of the area are arbitrary.
//   The additional working space in the suffix area is optional,
//   but the algorithm may be faster with a larger working space.
//   The maximum speed-up is achieved when the additional working
//   space equals the space for the suffixes.
//
// common_prefix_length is the length of a prefix that all the suffixes
//   to be sorted share (and which can be skipped in comparisons).
//
// target_prefix_length: A group of suffixes sharing a prefix of this
//   length do not have to be sorted further. On exit, such a group
//   is in the correct position with respect to other suffixes, but the
//   internal order of the group is arbitrary. If target_prefix_length
//   equals text length, the suffixes are completely sorted.
//   On the other hand, a smaller value keeps the worst case time
//   and space requirements reasonable (see below).
//
// report_finished_group is a functor that is called for every group
//   of suffixes that the algorithm is finished with (see
//   target_prefix_length). When StringsortSuffixes makes the call
//   report_finished_group(group_begin, group_end), it means that:
//   1. [group_begin,group_end) is a subrange of [suffix_area,suffixes_end)
//      and contains a subset of the original suffixes.
//   2. Either the size of the group is one or all suffixes in the group
//      share a prefix of length target_prefix_length.
//   3. The contents of [group_begin,group_end) are not changed or even
//      read anymore during the rest of the execution. Thus
//      report_finished_group may also change the contents.
//   Every suffix belongs to exactly one of the reported groups.
//   Use NullFinishedGroupReporter if no reporting is needed.
//
// TIME
//
// The time complexity of the algorithm is O(n log n + D)
// where n is the number of suffixes and D is the sum of the lengths
// of the distinguishing prefixes of the suffixes, i.e., the minimum
// number of characters that needs to be inspected to complete the job.
// In the worst case, D can be as large as n * (text length - (n-1)/2).
// Typically it is much smaller though.
// In particular, D <= n*target_prefix_length
//
// SPACE
//
// The algorithm uses little memory in addition to the input and the stack.
// However, the size of the stack can grow large. The size of the stack
// is at most a few hundred bytes times the following the quantity:
//   log(n) + min(target_prefix_length, F)
// F is the length of the longest prefix shared by more than
// kMinInsertionSortSize different suffixes.
// (kMinInsertionSortSize is defined in stringsort-inl.h)
////////////////////////////////////////////////////////////////
template <typename CharIterator, typename UnsignedInteger,
          typename FinishedGroupReporter>
UnsignedInteger* StringsortSuffixes(
    const CharIterator text, const CharIterator text_end,
    UnsignedInteger* suffix_area,
    UnsignedInteger* suffixes_in, UnsignedInteger* suffix_area_end,
    int64 common_prefix_length,
    int64 target_prefix_length,
    FinishedGroupReporter& report_finished_group);

// Use this as the FinishedGroupReporter if no reporting is needed.
// It is a bit faster than using your own null reporter.
template <typename UnsignedInteger>
class NullFinishedGroupReporter {
 public:
  void operator() (UnsignedInteger* begin, UnsignedInteger* end) const {}

  template <typename T>
  bool operator!=(const T&) const { return true; }
  template <typename T>
  bool operator!=(const NullFinishedGroupReporter<T>&) const { return false; }
};

}  // namespace dcsbwt

#endif  // DCSBWT_STRINGSORT_H__
