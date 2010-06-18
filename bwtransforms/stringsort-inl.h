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

// Implementation of suffix sorting using string sorting algorithms,
// i.e., algorithms that do not take advantage of the fact the input
// strings consists of suffixes of a single string.
//
// The main algorithm is string quicksort, a.k.a. multikey quicksort,
// with two different versions, one for the normal representation and
// one for the cached representation (see below). Small sets are
// sorted with insertion sort.
//
// The implementation consists of several functions that manipulate
// a single array containing the suffixes. On the initial call of
// StringsortSuffixes the array looks like this:
// ----------???????????????? (?=unsorted suffix -=empty space)
// On exit the array looks like this:
// ++++++++++++++++----------   (+=sorted suffix in its final place)
//
// During the execution the array looks like this:
// ++++++++----------????????
// or this:
// ++++++++------?#?#?#?#????   (#=cached characters)
// Each function call sees only a part of the array.
// On entry to the function that part looks like this:
//         ----------????
// or this:
//         ------?#?#?#?#
// When the function call exits, that part looks like this:
//         ++++----------
//
// Each function has essentially the same arguments as StringsortSuffixes
// and the explanations in stringsort.h are valid except for
// workarea, suffixes, and suffixes_end, which are pointers identifying
// the parts of the array on which the function operates, like this:
//         ----------????
//         ^         ^   ^
//         |         |   |
//   workarea  suffixes  suffixes_end
//
// Besides these three, only common_prefix_length changes during
// the execution. It is the length of a prefix shared by all the
// suffixes on which the function operates.
// The return value is a pointer to the start of the empty space:
//         ++++----------
//             ^
//             |
//             return value
//
// As indicated in the figures, there are two representations of the suffixes.
// The normal representation '?' is an integer, which is the starting position
// of the suffix in the text. The second representation '?#' consists of two
// integers. The first one is the starting position, and the second one
// is a cache containing characters of the suffix starting from the position
// common_prefix_length packed into the integer. A comparison of the caches
// is equivalent to comparing sizeof(UnsignedInteger) characters at once.
// The functions, whose names include Cache or Caching, operate on the
// representation with the cache and other functions operate on the normal
// representation.
//
// The purpose of caching is to reduce costly cache misses when accessing
// the text. The basic string quicksort algorithm accesses characters
// of a suffix one at a time and may take a long time between two
// accesses. This makes it likely that (the relevant part of) the suffix
// has been driven out of the system caches in the mean while.
// The caching algorithms, when they access the characters of a suffix,
// read multiple consecutive characters at a time and store them into
// the algorithmic cache next to the starting position.

#ifndef DCSBWT_STRINGSORT_INL_H__
#define DCSBWT_STRINGSORT_INL_H__

#include "stringsort.h"
#include "ternary_partition.h"

#include "inttypes.h"

#include <functional>  // for bind2nd, equal_to, greater_equal
#include <algorithm>   // for copy, copy_backward, find_if, partition
#include <cassert>

namespace dcsbwt {
namespace stringsort {

// Switch to insertion sort for groups no larger than this.
static const int kMaxInsertionSortSize = 20;

// Switch from normal representation to cached representation
// only for groups larger or equal to this.
static const int kMinCachingStartSize = 4096;

// Switch from cached representation to normal representation
// for groups smaller than this (but only at cache refill time).
static const int kMinCacheRefillSize = 400;

////////////////////////////////////////////////////////////////////
// Sort a set of suffixes (without caches) using a new
// character position, i.e., a new value of common_prefix_length.
// This is called for the initial input and every time
// common_prefix_length is incremented (when not caching).
// It does two special checks:
// 1. Has common_prefix_length reached target_prefix_length?
// 2. Is there a suffix of length common_prefix_length, i.e.,
//    a suffix for which the new position is past the end of the text?
////////////////////////////////////////////////////////////////////
template <typename CharIterator, typename UnsignedInteger,
          typename FinishedGroupReporter>
UnsignedInteger* StringsortNewPosition(
    const CharIterator text, const CharIterator text_end,
    UnsignedInteger* workarea,
    UnsignedInteger* suffixes, UnsignedInteger* suffixes_end,
    UnsignedInteger common_prefix_length,
    const UnsignedInteger target_prefix_length,
    FinishedGroupReporter& report_finished_group)
{
  const int64 num_suffixes = suffixes_end - suffixes;
  assert(num_suffixes > 0);
  // Optimization for the case of a single suffix
  if (1 == num_suffixes) {
    *workarea = *suffixes;
    report_finished_group(workarea, workarea + 1);
    return workarea + 1;
  }
  // Have we reached target_prefix_length?
  assert(common_prefix_length <= target_prefix_length);
  if (common_prefix_length == target_prefix_length) {
    if (workarea < suffixes) {
      suffixes_end = std::copy(suffixes, suffixes_end, workarea);
      suffixes = workarea;
    }
    report_finished_group(suffixes, suffixes_end);
    return suffixes_end;
  }
  // Handle the suffix of length common_prefix_length if it is there.
  UnsignedInteger finished_suffix = (text_end - text) - common_prefix_length;
  UnsignedInteger* finished_suffix_position
      = std::find(suffixes, suffixes_end, finished_suffix);
  if (finished_suffix_position != suffixes_end) {
    *finished_suffix_position = *suffixes;
    ++suffixes;
    *workarea = finished_suffix;
    report_finished_group(workarea, workarea + 1);
    ++workarea;
  }
  // Let StringsortDispatch do its magic for the remaining suffixes.
  return StringsortDispatch(text, text_end,
                            workarea, suffixes, suffixes_end,
                            common_prefix_length, target_prefix_length,
                            report_finished_group);
}

///////////////////////////////////////////////////////////////////
// Decide which algorithm to use for sorting.
// The candidates are:
// 1. InsertionSort for small sets.
// 2. Caching algorithms if conditions are right.
// 3. StringQuicksort otherwise.
///////////////////////////////////////////////////////////////////
template <typename CharIterator, typename UnsignedInteger,
          typename FinishedGroupReporter>
UnsignedInteger* StringsortDispatch(
    const CharIterator text, const CharIterator text_end,
    UnsignedInteger* workarea,
    UnsignedInteger* suffixes, UnsignedInteger* suffixes_end,
    UnsignedInteger common_prefix_length,
    const UnsignedInteger target_prefix_length,
    FinishedGroupReporter& report_finished_group)
{
  const int64 num_suffixes = suffixes_end - suffixes;
  assert(num_suffixes > 0);
  assert(common_prefix_length < target_prefix_length);
  // Check that we don't go past the end of the text for any suffix.
  assert(std::find_if(suffixes, suffixes_end,
                      std::bind2nd(std::greater_equal<UnsignedInteger>(),
                                   (text_end - text) - common_prefix_length))
         == suffixes_end);
  // Optimization for case of a single suffix
  if (1 == num_suffixes) {
    *workarea = *suffixes;
    report_finished_group(workarea, workarea + 1);
    return workarea + 1;
  }
  // Use insertion sort for small sets
  if (num_suffixes <= kMaxInsertionSortSize) {
    return InsertionSort(text, text_end,
                         workarea, suffixes, suffixes_end,
                         common_prefix_length, target_prefix_length,
                         report_finished_group);
  }
  // Use caching sort if the following conditions are met:
  // 1. There are enough (kMinCachingStartSize) suffixes to make
  //    caching worthwhile.
  // 2. There is enough workspace for the caches.
  // 3. We are still at least one cacheful away from reaching
  //    target_prefix_length.
  if (num_suffixes >= kMinCachingStartSize
     && suffixes - workarea >= num_suffixes
     && target_prefix_length - common_prefix_length >= sizeof(UnsignedInteger))
  {
    return StartCachingAndSort(text, text_end,
                               workarea, suffixes, suffixes_end,
                               common_prefix_length, target_prefix_length,
                               report_finished_group);
  }
  // Otherwise use the regular string quicksort
  return StringQuicksort(text, text_end,
                         workarea, suffixes, suffixes_end,
                         common_prefix_length, target_prefix_length,
                         report_finished_group);
}

///////////////////////////////////////////////////////////////////
// Functor for comparing two suffixes.
// Used by InsertionSort.
///////////////////////////////////////////////////////////////////
template <typename CharIterator, typename UnsignedInteger>
class SuffixLess : public std::binary_function<bool,
                                               UnsignedInteger,
                                               UnsignedInteger> {
 public:
  SuffixLess(CharIterator text, CharIterator text_end,
             UnsignedInteger common_prefix_length,
             UnsignedInteger target_prefix_length)
      : text_plus_prefix_(text+common_prefix_length), text_end_(text_end),
        distance_to_target_(target_prefix_length - common_prefix_length) {}
  bool operator() (UnsignedInteger a, UnsignedInteger b) const {
    CharIterator a_iter = text_plus_prefix_ + a;
    CharIterator b_iter = text_plus_prefix_ + b;
    assert(text_end_ - a_iter >= 0);
    assert(text_end_ - b_iter >= 0);
    const UnsignedInteger a_length = text_end_ - a_iter;
    const UnsignedInteger b_length = text_end_ - b_iter;
    // The comparison of a and b stops after limit characters
    // if no mismatch has been found by then.
    // This happens when the algorithm reaches the end of a,
    // the end of b or target_prefix_length, whichever comes first.
    UnsignedInteger limit
        = std::min<UnsignedInteger>(distance_to_target_, b_length);
    bool result_at_limit = false;
    if (a_length < limit) {
      limit = a_length;
      result_at_limit = true;
    }
    for (; limit; --limit) {
      unsigned char a_ch = *a_iter;
      unsigned char b_ch = *b_iter;
      if (a_ch != b_ch) return (a_ch < b_ch);
      ++a_iter; ++b_iter;
    }
    return result_at_limit;
  }
 private:
  const CharIterator text_plus_prefix_;
  const CharIterator text_end_;
  const UnsignedInteger distance_to_target_;
};

/////////////////////////////////////////////////////////////////
// Sort a (small) set of suffixes using standard insertion sort.
/////////////////////////////////////////////////////////////////
template <typename CharIterator, typename UnsignedInteger,
          typename FinishedGroupReporter>
UnsignedInteger* InsertionSort(
    const CharIterator text, const CharIterator text_end,
    UnsignedInteger* workarea,
    UnsignedInteger* suffixes, UnsignedInteger* suffixes_end,
    UnsignedInteger common_prefix_length,
    const UnsignedInteger target_prefix_length,
    FinishedGroupReporter& report_finished_group)
{
  const int64 num_suffixes = suffixes_end - suffixes;
  assert(num_suffixes > 0);
  // Optimization for case of a single suffix
  if (1 == num_suffixes) {
    *workarea = *suffixes;
    report_finished_group(workarea, workarea + 1);
    return workarea + 1;
  }
  // The actual insertion sort starts.
  // The suffixes will be moved from the end of the workarea
  // to the beginning of the workarea during sorting like this:
  // --------????   =>   ++------??   =>   ++++--------
  SuffixLess<CharIterator,UnsignedInteger> suffix_less(
      text, text_end,
      common_prefix_length, target_prefix_length);
  UnsignedInteger* begin = workarea;
  *begin = *suffixes++;
  UnsignedInteger* end = workarea + 1;
  while (suffixes != suffixes_end) {
    UnsignedInteger suffix = *suffixes++;
    UnsignedInteger* insertion_point = end;
    while (insertion_point > begin
           && suffix_less(suffix, *(insertion_point-1))) {
      *insertion_point = *(insertion_point-1);
      --insertion_point;
    }
    *insertion_point = suffix;
    ++end;
  }
  assert(end - begin == num_suffixes);
  // The insertion sort has finished.
  // Now [begin,end) contains the suffixes in order.
  // We still have to report the groups.
  if (NullFinishedGroupReporter<UnsignedInteger>() != report_finished_group) {
    for (UnsignedInteger* group_end = begin + 1;
         group_end < end; ++group_end) {
      if (suffix_less(*begin, *group_end)) {
        report_finished_group(begin, group_end);
        begin = group_end;
      }
    }
    report_finished_group(begin, end);
  }
  return end;
}

////////////////////////////////////////////////////////////////
// Functor for obtaining the character at position common_prefix_length
// for each suffix.
// Used by StringQuicksort.
////////////////////////////////////////////////////////////////
template <typename CharIterator, typename UnsignedInteger>
class SuffixChar : public std::unary_function<UnsignedInteger,unsigned char> {
 public:
  SuffixChar(CharIterator text, UnsignedInteger common_prefix_length)
      : text_plus_prefix_(text + common_prefix_length) {}
  unsigned char operator() (UnsignedInteger suffix) const {
    return text_plus_prefix_[suffix];
  }
 private:
  CharIterator text_plus_prefix_;
};

/////////////////////////////////////////////////////////////////
// Sort suffixes using string quicksort, a.k.a. multikey quicksort.
// Performs a ternary partition of the set using the characters
// at position common_prefix_length, the first position after the
// shared prefix. Recurses on each of the three parts.
/////////////////////////////////////////////////////////////////
template <typename CharIterator, typename UnsignedInteger,
          typename FinishedGroupReporter>
UnsignedInteger* StringQuicksort(
    const CharIterator text, const CharIterator text_end,
    UnsignedInteger* workarea,
    UnsignedInteger* suffixes, UnsignedInteger* suffixes_end,
    UnsignedInteger common_prefix_length,
    const UnsignedInteger target_prefix_length,
    FinishedGroupReporter& report_finished_group)
{
  assert(suffixes_end - suffixes > 0);
  SuffixChar<CharIterator,UnsignedInteger>
      suffix_key(text, common_prefix_length);
  UnsignedInteger* pivot = ChoosePivot(suffixes, suffixes_end, suffix_key);
  std::pair<UnsignedInteger*,UnsignedInteger*> equal_part
      = TernaryPartition(suffixes, suffixes_end, pivot, suffix_key);
  UnsignedInteger* equal_part_begin = equal_part.first;
  UnsignedInteger* equal_part_end = equal_part.second;
  if (suffixes < equal_part_begin) {
    workarea = StringsortDispatch(text, text_end,
                                  workarea, suffixes, equal_part_begin,
                                  common_prefix_length, target_prefix_length,
                                  report_finished_group);
  }
  assert(equal_part_end - equal_part_begin > 0);
  workarea = StringsortNewPosition(text, text_end,
                                   workarea, equal_part_begin, equal_part_end,
                                   common_prefix_length + 1,
                                   target_prefix_length,
                                   report_finished_group);
  if (equal_part_end < suffixes_end) {
    workarea = StringsortDispatch(text, text_end,
                                  workarea, equal_part_end, suffixes_end,
                                  common_prefix_length, target_prefix_length,
                                  report_finished_group);
  }
  return workarea;
}

////////////////////////////////////////////////////////////
// The representation of suffixes when caching is used.
////////////////////////////////////////////////////////////
template <typename UnsignedInteger>
struct CachedSuffix {
  UnsignedInteger suffix;
  UnsignedInteger cache;
  struct KeyGetter : public
        std::unary_function<CachedSuffix<UnsignedInteger>,UnsignedInteger> {
    UnsignedInteger operator() (
        CachedSuffix<UnsignedInteger> suffix) const {
      return suffix.cache;
    }
  };
  bool operator>= (const CachedSuffix<UnsignedInteger>& other) const {
    return this->suffix >= other.suffix;
  }
};

/////////////////////////////////////////////////////////////
// Fill the cache of all suffixes with new characters,
// i.e., with the characters in positions
// [common_prefix_length, common_prefix_length + sizeof(UnsignedInteger))
/////////////////////////////////////////////////////////////
template <typename CharIterator, typename UnsignedInteger>
void FillCache(
    const CharIterator text, const CharIterator text_end,
    CachedSuffix<UnsignedInteger>* suffixes,
    CachedSuffix<UnsignedInteger>* suffixes_end,
    const UnsignedInteger common_prefix_length)
{
  unsigned char buffer[sizeof(UnsignedInteger)];
  for (CachedSuffix<UnsignedInteger>* suffix = suffixes;
       suffix < suffixes_end; ++suffix) {
    // Copy characters to a local buffer.
    CharIterator begin = text + suffix->suffix + common_prefix_length;
    CharIterator end = begin + sizeof(UnsignedInteger);
    if (text_end < end) {
      end = text_end;
      memset(buffer, 0, sizeof(UnsignedInteger));
    }
    std::copy(begin, end, buffer);
    // Pack characters into a local integer in an endianness independent way.
    UnsignedInteger cache = buffer[0];
    for (int i = 1; i < sizeof(UnsignedInteger); ++i) {
      cache = (cache << 8) + buffer[i];
    }
    // Write cache to the array.
    suffix->cache = cache;
  }
}

////////////////////////////////////////////////////////////////
// Change from the standard representation (---------????)
// into the cached representation          (-----?#?#?#?#)
// and sort.
////////////////////////////////////////////////////////////////
template <typename CharIterator, typename UnsignedInteger,
          typename FinishedGroupReporter>
UnsignedInteger* StartCachingAndSort(
    const CharIterator text, const CharIterator text_end,
    UnsignedInteger* workarea,
    UnsignedInteger* uncached_suffixes, UnsignedInteger* uncached_suffixes_end,
    UnsignedInteger common_prefix_length,
    const UnsignedInteger target_prefix_length,
    FinishedGroupReporter& report_finished_group)
{
  const uint64 num_suffixes = uncached_suffixes_end - uncached_suffixes;
  assert(num_suffixes > kMaxInsertionSortSize);
  assert(uncached_suffixes - workarea >= num_suffixes);
  CachedSuffix<UnsignedInteger>* cached_suffixes_end
      = reinterpret_cast<CachedSuffix<UnsignedInteger>*>(
          uncached_suffixes_end);
  CachedSuffix<UnsignedInteger>* cached_suffixes
      = cached_suffixes_end - num_suffixes;
  // Copy starting positions to their new places.
  CachedSuffix<UnsignedInteger>* i = cached_suffixes;
  UnsignedInteger* j = uncached_suffixes;
  while (j < uncached_suffixes_end) {
    (i++)->suffix = *j++;
  }
  assert(i == cached_suffixes_end);
  // Fill the caches.
  FillCache(text, text_end,
            cached_suffixes, cached_suffixes_end, common_prefix_length);
  // Sort with string quicksort.
  return StringQuicksortWithCache(text, text_end,
                                  workarea,
                                  cached_suffixes, cached_suffixes_end,
                                  common_prefix_length, target_prefix_length,
                                  report_finished_group);
}

////////////////////////////////////////////////////////////////
// Change from the cached representation (?#?#?#?#)
// into the standard representation      (----????)
// Returns the range containing the uncached suffixes.
////////////////////////////////////////////////////////////////
template <typename UnsignedInteger>
std::pair<UnsignedInteger*,UnsignedInteger*> RemoveCache(
    CachedSuffix<UnsignedInteger>* cached_suffixes,
    CachedSuffix<UnsignedInteger>* cached_suffixes_end)
{
  assert(cached_suffixes_end - cached_suffixes > 0);
  UnsignedInteger* uncached_suffixes_end
      = reinterpret_cast<UnsignedInteger*>(cached_suffixes_end);
  UnsignedInteger* uncached_suffixes = uncached_suffixes_end;
  do {
    *--uncached_suffixes = (--cached_suffixes_end)->suffix;
  } while (cached_suffixes_end > cached_suffixes);
  return std::make_pair(uncached_suffixes, uncached_suffixes_end);
}

////////////////////////////////////////////////////////////////
// Refill caches if appropriate and sort.
// Before refilling, handle various special cases.
////////////////////////////////////////////////////////////////
template <typename CharIterator, typename UnsignedInteger,
          typename FinishedGroupReporter>
UnsignedInteger* RefillCacheAndSort(
    const CharIterator text, const CharIterator text_end,
    UnsignedInteger* workarea,
    CachedSuffix<UnsignedInteger>* suffixes,
    CachedSuffix<UnsignedInteger>* suffixes_end,
    UnsignedInteger common_prefix_length,
    const UnsignedInteger target_prefix_length,
    FinishedGroupReporter& report_finished_group)
{
  const uint64 num_suffixes = suffixes_end - suffixes;
  assert(num_suffixes > 0);
  assert(common_prefix_length <= target_prefix_length);
  // optimization for the case of one suffix
  if (num_suffixes == 1) {
    *workarea = suffixes->suffix;
    report_finished_group(workarea, workarea + 1);
    return workarea + 1;
  }
  // Find suffixes, whose end we have reached, i.e., for who the position
  // common_prefix_length is beyond the end of the text.
  UnsignedInteger first_finished_suffix
      = (text_end - text) - common_prefix_length;
  CachedSuffix<UnsignedInteger> finished_suffix;
  finished_suffix.suffix = first_finished_suffix;
  finished_suffix.cache = 0;
  CachedSuffix<UnsignedInteger>* cached_finished_end
      = std::partition(suffixes, suffixes_end,
             std::bind2nd(std::greater_equal<CachedSuffix<UnsignedInteger> >(),
                          finished_suffix));
  // If we found finished suffixes, handle them separately
  if (cached_finished_end > suffixes) {
    CachedSuffix<UnsignedInteger>* cached_finished_begin = suffixes;
    suffixes = cached_finished_end;
    // Now [cached_finished_begin,cached_finished_end) contains
    // the finished suffixes and [suffixes,suffixes_end) contains
    // the other suffixes.

    // Remove caches of the finished suffixes
    std::pair<UnsignedInteger*,UnsignedInteger*> uncached_finished
        = RemoveCache(cached_finished_begin, cached_finished_end);
    UnsignedInteger* finished_begin = uncached_finished.first;
    UnsignedInteger* finished_end = uncached_finished.second;
    // Sort in descending order of starting position.
    // This is their correct order.
    std::sort(finished_begin, finished_end,
              std::greater<UnsignedInteger>());
    // Move from the end of the workarea to the beginning of the workarea
    finished_end = std::copy(finished_begin, finished_end, workarea);
    finished_begin = workarea;
    workarea = finished_end;
    // Report the finished suffixes.
    for (UnsignedInteger* i = finished_begin; i != finished_end; ++i) {
      report_finished_group(i, i+1);
    }
    // If there were no unfinished suffixes, we are done.
    if (suffixes == suffixes_end) return workarea;
  }
  // If we have reached target_prefix_length, the  remaining suffixes
  // are finished, too. If so, they form a single group.
  if (common_prefix_length >= target_prefix_length) {
    std::pair<UnsignedInteger*,UnsignedInteger*> uncached
        = RemoveCache(suffixes, suffixes_end);
    UnsignedInteger* finished_suffixes = workarea;
    workarea = std::copy(uncached.first, uncached.second, workarea);
    report_finished_group(finished_suffixes, workarea);
    return workarea;
  }
  // Switch to non-caching algorithms if, either there are not
  // enough suffixes to justify caching, or there are less than
  // a cacheful of characters left before reaching target_prefix_length.
  if (num_suffixes < kMinCacheRefillSize
      || target_prefix_length - common_prefix_length < sizeof(UnsignedInteger))
  {
    std::pair<UnsignedInteger*,UnsignedInteger*> uncached_suffixes
        = RemoveCache(suffixes, suffixes_end);
    if (suffixes_end - suffixes <= kMaxInsertionSortSize) {
      return InsertionSort(text, text_end,
                           workarea,
                           uncached_suffixes.first, uncached_suffixes.second,
                           common_prefix_length, target_prefix_length,
                           report_finished_group);
    } else {
      return StringQuicksort(text, text_end,
                             workarea,
                             uncached_suffixes.first, uncached_suffixes.second,
                             common_prefix_length, target_prefix_length,
                             report_finished_group);
    }
  }
  // No more excuses left; refill the cache and sort.
  FillCache(text, text_end,
            suffixes, suffixes_end, common_prefix_length);
  return StringQuicksortWithCache(text, text_end,
                                  workarea, suffixes, suffixes_end,
                                  common_prefix_length, target_prefix_length,
                                  report_finished_group);
}

///////////////////////////////////////////////////////////////////
// Decide whether to continue cached sorting
// or to switch to non-caching insertion sort.
///////////////////////////////////////////////////////////////////
template <typename CharIterator, typename UnsignedInteger,
          typename FinishedGroupReporter>
UnsignedInteger* StringsortDispatchWithCache(
    const CharIterator text, const CharIterator text_end,
    UnsignedInteger* workarea,
    CachedSuffix<UnsignedInteger>* suffixes,
    CachedSuffix<UnsignedInteger>* suffixes_end,
    UnsignedInteger common_prefix_length,
    const UnsignedInteger target_prefix_length,
    FinishedGroupReporter& report_finished_group)
{
  const uint64 num_suffixes = suffixes_end - suffixes;
  assert(num_suffixes > 0);
  if (num_suffixes > kMaxInsertionSortSize) {
    return StringQuicksortWithCache(text, text_end,
                                    workarea, suffixes, suffixes_end,
                                    common_prefix_length, target_prefix_length,
                                    report_finished_group);
  }
  std::pair<UnsignedInteger*,UnsignedInteger*> uncached_suffixes
      = RemoveCache(suffixes, suffixes_end);
  return InsertionSort(text, text_end,
                       workarea,
                       uncached_suffixes.first, uncached_suffixes.second,
                       common_prefix_length, target_prefix_length,
                       report_finished_group);
}

///////////////////////////////////////////////////////////////
// String quicksort for the cached representation.
// This is the main cached sorting algorithm.
///////////////////////////////////////////////////////////////
template <typename CharIterator, typename UnsignedInteger,
          typename FinishedGroupReporter>
UnsignedInteger* StringQuicksortWithCache(
    CharIterator text, CharIterator text_end,
    UnsignedInteger* workarea,
    CachedSuffix<UnsignedInteger>* suffixes,
    CachedSuffix<UnsignedInteger>* suffixes_end,
    UnsignedInteger common_prefix_length, UnsignedInteger target_prefix_length,
    FinishedGroupReporter& report_finished_group)
{
  const uint64 num_suffixes = suffixes_end - suffixes;
  assert(num_suffixes > 0);
  typename CachedSuffix<UnsignedInteger>::KeyGetter get_key;
  CachedSuffix<UnsignedInteger>* pivot
      = ChoosePivot(suffixes, suffixes_end, get_key);
  std::pair<CachedSuffix<UnsignedInteger>*,CachedSuffix<UnsignedInteger>*>
      equal_part = TernaryPartition(suffixes, suffixes_end, pivot, get_key);
  CachedSuffix<UnsignedInteger>* equal_part_begin = equal_part.first;
  CachedSuffix<UnsignedInteger>* equal_part_end = equal_part.second;
  if (suffixes < equal_part_begin) {
    workarea = StringsortDispatchWithCache(text, text_end,
                                   workarea, suffixes, equal_part_begin,
                                   common_prefix_length, target_prefix_length,
                                   report_finished_group);
  }
  assert(equal_part_end - equal_part_begin > 0);
  workarea = RefillCacheAndSort(text, text_end,
                                workarea, equal_part_begin, equal_part_end,
                                common_prefix_length + sizeof(UnsignedInteger),
                                target_prefix_length,
                                report_finished_group);
  if (equal_part_end < suffixes_end) {
    workarea = StringsortDispatchWithCache(text, text_end,
                                   workarea, equal_part_end, suffixes_end,
                                   common_prefix_length, target_prefix_length,
                                   report_finished_group);
  }
  return workarea;
}

}  // namespace stringsort

/////////////////////////////////////////////////////////////////////
// The interface function, see stringsort.h
/////////////////////////////////////////////////////////////////////
template <typename CharIterator, typename UnsignedInteger,
          typename FinishedGroupReporter>
UnsignedInteger* StringsortSuffixes(
    const CharIterator text, const CharIterator text_end,
    UnsignedInteger* suffix_area, UnsignedInteger* suffixes_in,
    UnsignedInteger* suffix_area_end,
    int64 common_prefix_length,
    int64 target_prefix_length,
    FinishedGroupReporter& report_finished_group)
{
  assert(text_end - text >= 0);
  assert(suffixes_in - suffix_area >= 0);
  assert(suffix_area_end - suffixes_in >= 0);
  assert(common_prefix_length >= 0);
  assert(target_prefix_length >= common_prefix_length);
  // Check that target_prefix_length (and consequently common_prefix_length)
  // is a value that can be represented by UnsignedInteger.
  assert(static_cast<int64>(static_cast<UnsignedInteger>(
                             target_prefix_length)) >= target_prefix_length);
  if (suffixes_in == suffix_area_end) {
    return suffix_area;
  }
  UnsignedInteger* suffixes_end =
      stringsort::StringsortNewPosition(
          text, text_end,
          suffix_area, suffixes_in, suffix_area_end,
          static_cast<UnsignedInteger>(common_prefix_length),
          static_cast<UnsignedInteger>(target_prefix_length),
          report_finished_group);
  assert(suffixes_end - suffix_area == suffix_area_end - suffixes_in);
  return suffixes_end;
}

}  // namespace dcsbwt

#endif  // DCSBWT_STRINGSORT_INL_H__
