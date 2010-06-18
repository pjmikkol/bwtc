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

/***************************************************************************
 * A difference cover for a period length p is a set D of integers         *
 * that satisfies: { (i-j) mod p | i,j in D } == [0,p)                     *
 * In other words, the differences of the integers (mod p)                 *
 * cover the whole period.                                                 *
 *                                                                         *
 * Example: D=={0,1,3,7,8} is a difference cover for period length 16      *
 *   0-0 == 0    8-0       ==  8                                           *
 *   1-0 == 1    0-7 == -7 ==  9 (mod 16)                                  *
 *   3-1 == 2    1-7 == -6 == 10 (mod 16)                                  *
 *   3-0 == 3    3-8 == -5 == 11 (mod 16)                                  *
 *   7-3 == 4    3-7 == -4 == 12 (mod 16)                                  *
 *   8-3 == 5    0-3 == -3 == 13 (mod 16)                                  *
 *   7-1 == 6    1-3 == -2 == 14 (mod 16)                                  *
 *   7-0 == 7    0-1 == -1 == 15 (mod 16)                                  *
 *                                                                         *
 * This file defines two classes to support the use of difference          *
 * covers in the suffix array construction and Burrows-Wheeler transform   *
 * algorithms described in the following articles (available at            *
 * http://www.cs.helsinki.fi/juha.karkkainen/):                            *
 *                                                                         *
 *   J. Karkkainen, S. Burkhardt:                                          *
 *   Fast lightweight suffix array construction and checking.              *
 *   In Proc. 14th Symposium on Combinatorial Pattern Matching (CPM '03).  *
 *   LNCS 2676, Springer, 2003, pp. 55-69.                                 *
 *                                                                         *
 *   J Karkkainen, P. Sanders, S. Burkhardt:                               *
 *   Linear work suffix array construction.                                *
 *   J. ACM, 53 (6), pp. 918-936, 2006.                                    *
 *                                                                         *
 *   J. Karkkainen:                                                        *
 *   Fast BWT in Small Space by Blockwise Suffix Sorting.                  *
 *   Theoretical Computer Science, 2007 (to appear).                       *
 *                                                                         *
 * The class DifferenceCover, in particuler, is quite general and          *
 * may be useful for other applications of difference covers.              *
 ***************************************************************************/

#ifndef DCSBWT_DIFFERENCE_COVER_H__
#define DCSBWT_DIFFERENCE_COVER_H__

#include <vector>
#include <cassert>

namespace dcsbwt {

//-------------------------
class DifferenceCover {
  // The class DifferenceCover computes and stores a difference cover D
  // for a period length p, where p is any power of two representable by
  // unsigned int. The size of D for small values of p is:
  //   p   1  2  4  8  16  32  64  128  256  512 1024 2048 4096 8192 16386
  //  |D|  1  2  3  4   5   7   9   13   20   28   40   58   82  112   160
  // For larger p, the size of D is guaranteed to be smaller than
  // sqrt(1.5*p)+5.5. (The size of the supporting data structures is
  // more than 8*p bytes, though, so be careful with large values of p).
  //
  // Note: unsigned int is used in most computations, because we want
  // to do arithmetic modulo a power of two (see C++ Standard 3.9.1/4).
  // It is possible to use another integral type, even a user-defined type
  // that behaves like an integer, for providing the input for and storing
  // the output of the member functions as long as the following requirements
  // are satisfied:
  // - The type can represent all values in the range [0,period]
  //   (note: including period) and conversion between the type and
  //   unsigned int works both ways for these values.
  // - For any value, a conversion from the type to unsigned int is
  //   performed mod 2^w, where w is the number of bits in unsigned int
  //   (see C++ Standard 4.7/2). For example, -1 is converted to 2^w-1.
  //   (This concerns only ModuloPeriod.)
  // - The type supports division by a power of two using operator>>
  //   (for any value used as an input to DivideByPeriod).
  // All built-in integral types satisfy the requirements.

 public:
  // constructor            require: period is a power of two
  explicit DifferenceCover(unsigned int period);
  // no default constructor
  // implicit destructor
  // implicit copy constructor
  // implicit copy assignment

  // get basic dimensions
  unsigned int period() const { return 1 << logperiod_; }
  unsigned int size() const { return cover_.size(); }

  // Division by period: i div p and i mod p
  // The results satisfy: (i div p) * p + (i mod p) == i
  // See comments above on the type of i (for both div and mod).
  // The type Integer must support division by a power of two using
  // operator>> for the value i. For built-in types, the behaviour
  // for negative values is implementation defined (5.8/3).
  // Thus negative inputs are forbidden.
  template <typename Integer>
  Integer DivideByPeriod (Integer i) const {
    assert(i >= 0);
    return i >> logperiod_;
  }
  // ModuloPeriod accepts any value as an argument, including negative values
  // and values larger than MAX_UINT, provided the conversion is done
  // mod 2^w (see above). The return value is always in the range [0,period).
  // For example, ModuloPeriod(-1) returns period-1.
  unsigned int ModuloPeriod (unsigned int i) const { return i & mask_; }

  // returns true if i is in D
  // require: i < p
  inline bool Contains(unsigned int i) const;

  // rank/select
  // Rank(i) = size of (D intersect [0,i))               require: i<p
  // Select(j) = jth element of D (starting with 0th)    require: j<|D|
  // Invariant: Rank(Select(j))==j
  //            Select(Rank(i))>=i when Rank(i)<|D|
  inline unsigned int Rank(unsigned int i) const;
  inline unsigned int Select(unsigned int j) const;

  // returns k in D such that (k+diff) mod p is in D too
  // require: diff < p
  inline unsigned int Coverer(unsigned int diff) const;

  // STL random access iterators for D
  typedef std::vector<unsigned int>::const_iterator iterator;
  iterator begin() const { return cover_.begin(); }
  iterator end() const { return cover_.end(); }

  // Find the size of a difference cover without constructing one.
  static unsigned int CoverSizeForPeriod(unsigned int period);

  // Helps keep track of space requirements.
  static size_t SizeInBytesForPeriod(unsigned int period);

 private:
  const unsigned int logperiod_;
  const unsigned int mask_;         // period-1 == 00..011..1

  std::vector<unsigned int> cover_;
  std::vector<unsigned int> ranks_;
  std::vector<unsigned int> coverers_;

  // Most of the construction is done in these functions.
  void ComputeCover();
  void ComputeRanks();
  void ComputeCoverers();

  // white box test
  bool IsCorrect() const;
};

template <typename Integer>
class DifferenceCoverSample {
  // The class DifferenceCoverSample represents the set
  // Dn = { 0<=i<n | i mod p is in D }
  // where D is a difference cover for a period p.
  // The key property of Dn is the support for the following operation:
  // - Shift(i,j) returns k in [0,p) such that i+k and j+k
  //   are both in Dn (for any i,j in [0,n-k))
  // Other operations provided by the class include:
  // - Fill(OutputIterator) writes Dn to the output range
  // - Rank(i) (for i in Dn) is the position of i in the Fill output
  // - PeriodInterval() is the distance from i to i+p in the Fill output
  //   (must be the same for all i, see note on ordering below)
  //
  // Note: The ordering of Dn produced by Fill must satisfy the
  // condition that Rank(i+p)-Rank(i) has the same (positive)
  // value for all i in Dn with i<n-p. The function PeriodInterval
  // returns this value. Two orders satisfying this are
  // (d0<d1<d2<... are the elements of D):
  // 1: d0, d1, d2, ..., d0+p, d1+p, d2+p, ..., d0+2*p, ...
  //   This increasing order  is simple to implement
  //   and may have an advantage in the speed of Rank.
  // 2: d0, d0+p, d0+2*p, ..., d1, d1+p, d1+2*p, ..., d2, ...
  //   This order corresponds to the one in the articles mentioned above
  //   and may make sorting the sample suffixes simpler and faster.
  // The current implementation uses the first (increasing) order.
  //
  // The template argument Integer of can be any integral type, even
  // a user defined type that behaves like an integer as long as
  // the following requirements are satisfied:
  // - conversion to unsigned int should be performed mod 2^w
  //   where w is the number of bits in unsigned int (see 4.7/2)
  // - support the following operations when all values and results
  //   are in the range [0,n]:
  //   + conversion from unsigned int
  //   + division and multiplication by a power of two using
  //     operator>> and operator<<
  //   + addition (+), subtraction (-), multiplication (*), pre-increment (++)
  //   + order comparison (<)
  // All built-in integral types satisfy the requirements

 public:
  typedef Integer integer_type;

  // constructor
  // require: period is a power of two
  // require: range >= period
  DifferenceCoverSample(unsigned int period, Integer range);
  // no default constructor
  // implicit destructor
  // implicit copy constructor
  // implicit copy assignment

  // get basic dimensions
  inline Integer period() const { return period_; }          // the p
  inline Integer range() const { return range_; }            // the n
  inline Integer size() const { return size_; }              // size of Dn
  inline Integer period_size() const { return dc_.size(); }  // size of D

  // Contains(i) returns true if i is in Dn
  inline bool Contains(Integer i) const;

  // Shift (see above)
  // The type of arguments is unsigned int since only the values
  // i mod p and j mod p are used in the computation.
  // Note: i and j can be any values. No check is made to ensure
  // that i, j, i+shift, or j+shift is in the range.
  inline Integer Shift(unsigned int i, unsigned int j) const;

  // Fill an output range with the elements of Dn.
  // The output range must have room for size() elements.
  // The value type of the output range must be able to store
  // integers in the range [0,range()).
  template <typename OutputIterator>
  OutputIterator Fill(OutputIterator it) const;

  // Rank (see above)
  // require: i is in Dn
  inline Integer Rank(Integer i) const;

  // PeriodInterval (see above)
  inline Integer PeriodInterval() const;

 private:

  const DifferenceCover dc_;
  const Integer period_;
  const Integer range_;
  Integer size_;

  void ComputeSize();
};

}  // namespace dcsbwt

#endif  // DCSBWT_DIFFERENCE_COVER_H__
