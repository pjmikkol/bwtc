/***********************************************************************
 * Implementation is originally from:                                  *
 * http://code.google.com/p/dcs-bwt-compressor/                        *
 ***********************************************************************/

#include "difference_cover-inl.h"
#include "difference_cover.h"

#include <cassert>     // for assert
#include <iterator>    // for back_insert_iterator
#include <algorithm>   // for copy, fill_n, adjacent_find, find_if
#include <functional>  // for bind2nd, greater_equal, not_equal_to
#include <numeric>     // for partial_sum

namespace dcsbwt {

namespace {

// Precomputed difference covers
// They were computed using a branch-and-bound search.
// All are minimal difference covers for these periods except the one
// for period 256, for which the search was not completed.

unsigned int const cover1[] = {0};
unsigned int const cover2[] = {0,1};
unsigned int const cover4[] = {0,1,3};
unsigned int const cover8[] = {0,1,3,4};
unsigned int const cover16[] = {0,1,3,7,8};
unsigned int const cover32[] = {0,1,3,7,12,17,25};
unsigned int const cover64[] = {0,1,3,8,18,34,40,44,53};
unsigned int const cover128[] = {0,1,3,7,17,40,55,64,75,85,104,109,117};
unsigned int const cover256[] = {
  0,1,3,7,12,20,30,44,65,80,89,96,114,122,128,150,196,197,201,219};
unsigned int const * const precomputed_covers[] = {
  cover1,
  cover2,
  cover4,
  cover8,
  cover16,
  cover32,
  cover64,
  cover128,
  cover256
};
unsigned int const precomputed_sizes[] = {
  1, 2, 3, 4, 5, 7, 9, 13, 20
};


// Base two logarithm of a power of two.
// Fails if the input is not a power of two.
unsigned int Log2(unsigned int power_of_two) {
  unsigned int log = 0;
  unsigned int p = (power_of_two >> 1);
  while (p) {
    ++log;
    p >>= 1;
  }
  // fail if the input was not a power of two
  assert((1U << log) == power_of_two);
  return log;
}

// find if there is a precomputed cover
bool ExistsPrecomputedCover(unsigned int period) {
  return Log2(period) <= Log2(256);
}

// compute a parameter for the Colbourn-Ling algorithm (see below
unsigned int ColbournLingDegree(unsigned int period) {
  unsigned int r = 0;
  while (24*r*r+36*r+13 < period) ++r;
  return r;
}

} // unnamed namespace

//----------------------------------
// DifferenceCover member functions
//----------------------------------

// static member function for computing the cover size without construction
unsigned int DifferenceCover::CoverSizeForPeriod(unsigned int period) {
  if (ExistsPrecomputedCover(period)) {
    return precomputed_sizes[Log2(period)];
  } else {
    return 4 + 6 * ColbournLingDegree(period);
  }
}

// static member function for computing the object size without construction
size_t DifferenceCover::SizeInBytesForPeriod(unsigned int period) {
  return 3 * sizeof(std::vector<unsigned int>) +
      (2 * static_cast<size_t>(period) +
       CoverSizeForPeriod(period) + 2) * sizeof(unsigned int);
}

// constructor
DifferenceCover::DifferenceCover(unsigned int period)
    : logperiod_(Log2(period)), mask_(period-1) /*TODO!*/ {
  ComputeCover();
  ComputeRanks();
  ComputeCoverers();
  IsCorrect();
}

// compute the difference cover
// for small periods use precomputed cover
// for large periods use the algorithm described in
//    C.J. Colbourn and A.C.H. Ling:
//    Quorums from Difference Covers.
//    Information Processing Letters 75(1-2):9-12, 2000.
void DifferenceCover::ComputeCover() {
  if (ExistsPrecomputedCover(period())) {
    // use precomputed cover
    unsigned int size = precomputed_sizes[logperiod_];
    cover_.reserve(size);
    unsigned int const * cp = precomputed_covers[logperiod_];
    std::back_insert_iterator<std::vector<unsigned int> > it(cover_);
    std::copy(cp, cp+size, it);
  } else {
    // use the Colbourn-Ling algorithm
    unsigned int r = ColbournLingDegree(period());
    unsigned int size = 4 + 6*r;
    cover_.reserve(size);
    std::back_insert_iterator<std::vector<unsigned int> > it(cover_);
    *it++ = 0;
    std::fill_n(it, r, 1);
    *it++ = r+1;
    std::fill_n(it, r, 2*r+1);
    std::fill_n(it, 2*r+1, 4*r+3);
    std::fill_n(it, r+1, 2*r+2);
    std::fill_n(it, r, 1);
    std::partial_sum(cover_.begin(), cover_.end(), cover_.begin());

    // In general, at this point the largest (last) elements could
    // be larger than the period length, and then one would have to
    // take modulo period, re-sort the elements and remove duplicates.
    // However, if period > 90 then all elements are in the period.
    // We start using Colbourn-Ling from period length 256, so just
    // do a check here to be safe.
    assert(cover_.back() < period());
  }
}

// Compute the lookup table for finding ranks.
// (See IsCorrect below.)
void DifferenceCover::ComputeRanks() {
  ranks_.reserve(period());
  std::back_insert_iterator<std::vector<unsigned int> > inserter(ranks_);
  std::vector<unsigned int>::const_iterator it;
  *inserter++ = 0;
  unsigned int pos = 0;
  for (it = cover_.begin(); it != cover_.end(); ++it) {
    std::fill_n(inserter, (*it)-pos, it-cover_.begin());
    pos = *it;
  }
  std::fill_n(inserter, period()-pos-1, size());
  assert(ranks_.size() == period());
}

// Compute lookup table for finding cover elements for a given difference.
// (See IsCorrect below.)
void DifferenceCover::ComputeCoverers() {
  coverers_.resize(period());
  coverers_[0] = cover_.front();
  std::vector<unsigned int>::const_iterator it1, it2;
  for (it1 = cover_.begin(); it1 != cover_.end(); ++it1) {
    for (it2 = cover_.begin(); it2 != it1; ++it2) {
      assert(*it1 > *it2);
      unsigned int diff = (*it1 - *it2);
      coverers_[diff] = *it2;
      coverers_[period()-diff] = *it1;
    }
  }
}

// Verify the correctness of the computed tables.
bool DifferenceCover::IsCorrect() const {

  // cover_ should contain distinct integers from [0,period)
  // in increasing order.
  // Check increasing order.
  assert(cover_.end() ==
           std::adjacent_find(cover_.begin(), cover_.end(),
                              std::greater_equal<unsigned int>() ));
  // Check that they are in [0,period)
  assert(cover_.back() < period());

  // ranks_ should contain 0s up to and including the first cover element,
  // 1s up to and including the second element, etc.
  // For example, if period==16:
  //   cover_ 0 1   3       7 8
  //   ranks_ 0 1 2 2 3 3 3 3 4 5 5 5 5 5 5 5
  assert(ranks_.size() == period());
  std::vector<unsigned int>::const_iterator it = ranks_.begin();
  for (unsigned int i = 0; i < size(); ++i) {
    it = std::find_if(it, ranks_.end(),
                      std::bind2nd(std::not_equal_to<unsigned int>(), i) );
    assert(static_cast<unsigned int>(it-ranks_.begin()) == cover_[i]+1);
  }
  it = std::find_if(it, ranks_.end(),
                    std::bind2nd(std::not_equal_to<unsigned int>(), size()) );
  assert(ranks_.end() == it);

  // coverers_ should contain at each position diff in [0,period)
  // an integer i such that both i and (i+diff) mod period are in
  // the cover.
  // For example, if period==16 and _cover=={0,1,3,7,8}
  //                 diff  0 1 2 3 4 5 6 7 8  9 10 11 12 13 14 15
  //      _coverers[diff]  0 0 1 0 3 3 1 0 0  7  7  8  7  3  3  1
  // _coverers[diff]+diff  0 1 3 3 7 8 7 7 8 16 17 19 19 16 17 16
  // (_coverers[diff]+diff) mod period        0  1  3  3  0  1  0
  assert(coverers_.size() == period());
  for (unsigned int diff = 0; diff < period(); ++diff) {
    unsigned int i = coverers_[diff];
    assert(i < period());
    unsigned int j = ModuloPeriod(i + diff);
    assert(Contains(i));
    assert(Contains(j));
  }

  return true;
}

}  // namespace dcsbwt
