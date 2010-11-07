/**
 * @file sequence_detector.cc
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
 * Implementation of SequenceDetector-class. Detection of repeating
 * sequencies.
 */

#include <cassert>
#include "sequence_detector.h"
#include "preprocessor.h"
#include "../globaldefs.h"

namespace bwtc {
namespace long_sequences {

namespace prime {

int mulm(int a, int b, int m)
{
  return static_cast<int64>(a)*b % m;
}

} // namespace prime

bool IsPrime(int n)
{
  using namespace prime;
  assert(n % 2 == 1);
  if ( n < 1373653) return (!Witness(2,n) && !Witness(3,n));
  else if ( n < 9080191 ) return (!Witness(31,n) && !Witness(73,n));
  /* n < 4,759,123,141 */
  else return (!Witness(2,n) && !Witness(7,n) && !Witness(61,n));
}

int FindPrimeLeq(int n) {
  if(n <= 3) return 3;
  assert(n > 3);
  if ((n&1) == 0) n -= 1;
  while (!IsPrime(n)) n -= 2;
  return n;
}

namespace hash_functions {

PrimeHasher::PrimeHasher() : prev_hash_(0) {}

unsigned PrimeHasher::Initialize(uint32 size_of_table, uint32 window_length)
{
  q_ = FindPrimeLeq(size_of_table);
  c_ = 1;
  for(unsigned i = 0; i < window_length - 1; ++i) {
    c_ <<= 8;
    c_ %= q_;
  }
  return q_;
}

int64 PrimeHasher::InitValue(const byte *from, unsigned len) {
  prev_hash_ = 0;
  for(unsigned i = 0; i < len; ++i)
    prev_hash_ = ((prev_hash_ << 8) + *from++) % q_;
  return prev_hash_;
}

int64 PrimeHasher::Update(byte old_val, byte new_val) {
  prev_hash_ = (((prev_hash_ - old_val*c_) << 8) + new_val) % q_;
  if (prev_hash_ < 0) prev_hash_ += q_;
  return prev_hash_;
}

/* Used in hashing */
const uint32 MAGIC = 37;
//const uint32 MAGIC = 257;

MaskHasher::MaskHasher() : prev_hash_(0) {}

unsigned MaskHasher::Initialize(uint32 size_of_table, uint32 window_length)
{
  //TODO: Put more care for choosing the mask_ and size of hash table
  mask_ = MostSignificantBit(size_of_table) - 1;
  c_ = 1;
  for(unsigned i = 1; i < window_length; ++i)
    c_ = (c_*MAGIC) & mask_;
  return mask_ + 1;
}

int64 MaskHasher::InitValue(const byte *from, unsigned len) {
  prev_hash_ = *from++;
  for(unsigned i = 1; i < len; ++i)
    prev_hash_ = ((prev_hash_*MAGIC) + *from++) & mask_;
  return prev_hash_;
}

int64 MaskHasher::Update(byte old_val, byte new_val) {
  prev_hash_ = ((prev_hash_ - old_val*c_)*MAGIC + new_val) & mask_;
  return prev_hash_;
}
} // namespace hash_functions


} // namespace long_sequences
} // namespace bwtc
