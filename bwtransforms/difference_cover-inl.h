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

#ifndef DCSBWT_DIFFERENCE_COVER_INL_H__
#define DCSBWT_DIFFERENCE_COVER_INL_H__

#include "difference_cover.h"

#include <cassert>     // for assert

#include <algorithm>   // for transform
#include <functional>  // for bind2nd, plus

namespace dcsbwt {

//------------------
// DifferenceCover
//------------------

inline bool DifferenceCover::Contains(unsigned int i) const {
  assert(i < period());
  return Rank(i) < size() && Select(Rank(i)) == i;
}
inline unsigned int DifferenceCover::Rank(unsigned int i) const {
  assert(i < period());
  return ranks_[i];
}
inline unsigned int DifferenceCover::Select(unsigned int j) const {
  assert(j < size());
  return cover_[j];
}
inline unsigned int DifferenceCover::Coverer(unsigned int diff) const {
  assert(diff < period());
  return coverers_[diff];
}

//-------------------------
// DifferenceCoverSample
//-------------------------

// constructor
template <typename Integer>
DifferenceCoverSample<Integer>::DifferenceCoverSample(
    unsigned int period, Integer range )
    : dc_(period), period_(period), range_(range)
{
  // Check that range >= period
  // We cannot convert range into unsigned int since its value
  // could legally exceed UINT_MAX. It is also possible that
  // the value of period exceeds the maximum value of Integer,
  // but then period > max(Integer) >= range, which violates
  // the condition. Thus one of the following tests fails if
  // and only if period > range.
  assert(static_cast<unsigned int>(static_cast<Integer>(period)) == period);
  assert(range >= static_cast<Integer>(period));
  ComputeSize();
}

template <typename Integer>
void DifferenceCoverSample<Integer>::ComputeSize() {
  Integer n_div_p = dc_.DivideByPeriod(range_);
  unsigned int n_mod_p = dc_.ModuloPeriod(range_);
  size_ = static_cast<Integer>(dc_.size()) * n_div_p
      + static_cast<Integer>(dc_.Rank(n_mod_p));
}

template <typename Integer>
inline bool DifferenceCoverSample<Integer>::Contains(Integer i) const {
  return 0 <= i && i < range() && dc_.Contains(dc_.ModuloPeriod(i));
}

template <typename Integer>
inline Integer DifferenceCoverSample<Integer>::Shift(
    unsigned int i, unsigned int j) const
{
  unsigned int diff_mod_p = dc_.ModuloPeriod(j-i);
  return dc_.ModuloPeriod(dc_.Coverer(diff_mod_p) - i);
}

// Fill the output range with the elements of Dn
// Dn order 1 (see header)
template <typename Integer>
template <typename OutputIterator>
OutputIterator DifferenceCoverSample<Integer>::Fill(
    OutputIterator out_iter) const
{
  Integer start_of_period = 0;
  Integer num_elements_written = 0;
  // Write copies of the difference cover, each shifted by one
  // period length from the previous one.
  // Note that the test in while is equivalent to
  // start_of_period + period() <= range(), but the latter could
  // overflow, since Integer is not guaranteed to be able to represent
  // a value larger than range().
  while (start_of_period <= range() - period()) {
    out_iter = std::transform(dc_.begin(), dc_.end(), out_iter,
                   std::bind2nd(std::plus<Integer>(), start_of_period));
    start_of_period += period();
    num_elements_written += dc_.size();
  }
  // Now there might not be room for a full copy of the difference
  // cover anymore, so write the last copy one element at a time.
  DifferenceCover::iterator dc_iter = dc_.begin();
  while ((dc_iter != dc_.end()) &&
         (start_of_period < range() - static_cast<Integer>(*dc_iter)) ) {
    *out_iter++ = start_of_period + *dc_iter++;
    ++num_elements_written;
  }
  assert(size() == num_elements_written);
  return out_iter;
}

// Rank
// Dn order 1 (see header)
template <typename Integer>
inline Integer DifferenceCoverSample<Integer>::Rank(Integer i) const {
  assert(Contains(i));
  Integer i_div_p = dc_.DivideByPeriod(i);
  unsigned int i_mod_p = dc_.ModuloPeriod(i);
  Integer num_elements_in_earlier_periods = period_size() * i_div_p;
  Integer num_earlier_elements_in_this_period = dc_.Rank(i_mod_p);
  return num_elements_in_earlier_periods + num_earlier_elements_in_this_period;
}

// period_distance
// Dn order 1 (see header)
template <typename Integer>
inline Integer DifferenceCoverSample<Integer>::PeriodInterval() const {
  return dc_.size();
}

}  // namespace dcsbwt

#endif  // DCSBWT_DIFFERENCE_COVER_INL_H__
