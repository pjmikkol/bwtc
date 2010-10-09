/**
 * @file sequence_detector.h
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
 * Header for structs and functions that are implemented and used in
 * SequenceDetector-class.
 */

#ifndef BWTC_SEQUENCE_DETECTOR_H_
#define BWTC_SEQUENCE_DETECTOR_H_

#include <list>
#include "../globaldefs.h"

namespace bwtc {
namespace long_sequences {

/**
 * Utilities for testing prime numbers. Test we use here is deterministic
 * variant of Miller-Rabin primality test.
 */
namespace prime {

/**
 * Product of a and b in modulo m
 * @param a member of multiplication
 * @param b other member of multiplication
 * @param m modulo where the multiplication is conducted
 * @return (a*b) % m
 */
int mulm(int a, int b, int m);

/**
 * Computes b to the power of e in modulo m.
 * @param b base
 * @param e exponent
 * @param m modulo where the exponentiation is conducted.
 * @return (b**e) % m
 */
template<class T>
T powm(T b, T e, T m)
{
  if (e == 0) return 1;
  T x = powm(b, e>>1, m);
  T xx = mulm(x, x, m);
  return (e & 1) ? mulm(b, xx, m) : xx;
}

/**
 * Checks if second argument is composite with Miller-Rabin primality test.
 * There is small probability for mistakes.
 *
 * @param a possible witness fo compositeness of n
 * @param n number to be inspected
 * @return true if n is composite for sure
 */
template<class T>
bool Witness(T a, T n)
{
  // TODO: Replace use of builtin for other compilers than gcc
  int t = __builtin_ffsll(n-1)-1;
  T x = powm(a, (n-1) >> t, n);
  for(int i = 0; i < t; ++i) {
    T y = mulm(x,x,n);
    if (y == 1 && x != 1 && x != n-1) return 1;
    x = y;
  }
  return x != 1;
}
} //namespace prime

/**
 * Determinized version of Miller-Rabin primality test.
 *
 * @param n number to be inspected
 * @return true if n is prime, false if composite
 */
bool IsPrime(int n);

/**
 * Finds prime which is lesser or equal to given n.
 *
 * @param n upper bound for prime returned
 * @return prime lesser or equal to n
 */
int FindPrimeLeq(int n);

/**
 * Classes for computing rolling hash functions.
 *
 * Each Hasher-object is designed to used as a template parameter and is
 * thought to be inherited from idealized Hasher-class. We don't use
 * class-hierarchy since we want to avoid overhead from virtual function calls.
 *
 * Hasher-classes need to implement the following functions:
 * unsigned Initialize(unsigned size_of_table);
 * int64 InitValue(const byte *from, unsigned length);
 * int64 Update(byte old_val, byte new_val);
 */
namespace hash_functions {

/**
 * Computes rolling hash using remainder of division with primes.
 */
class PrimeHasher {
 public:
  PrimeHasher();
  /**
   * Initializes the hasher and seeks parameters used for calculating hash
   * values. Client must supply constructor with the maximum size of the
   * hash table.
   *
   * @param size_of_table maximum size of hash table
   * @param window_length length of the rolling hash
   * @return size of hash table
   */
  unsigned Initialize(unsigned size_of_table, uint32 window_length);
  /**
   * Initializes the value of rolling hash with the hash value of
   * [*from, *(from +len) ).
   *
   * @param from pointer to the start of sequence
   * @param len length of sequence
   * @return hash value of sequence [*from, *(from +len) )
   */
  int64 InitValue(const byte *from, unsigned len);
  /**
   * Get the expected size of the hash table.
   *
   * @return expected size of the hash table
   */
  unsigned Size() {
    return q_;
  }
  /**
   * Updates the hash value from the hash of [old_val, new_val) to
   * (old_val, new_val].
   *
   * @return updated hash value
   */
  int64 Update(byte old_val, byte new_val);
 private:
  /**<Hash value of previous sequence*/
  int64 prev_hash_; 
  /**<Parameter used in calculation of hash function*/
  int64 c_; 
  /**<Prime used in the calculation of hash function. Hash values will be in
     range [0,q_) */
  int64 q_; 
};


/**
 * Computes rolling hash using bitmask.
 */
class MaskHasher {
 public:
  MaskHasher();
  /**
   * Initializes the hasher and seeks parameters used for calculating hash
   * values. Client must supply constructor with the maximum size of the
   * hash table.
   *
   * @param size_of_table maximum size of hash table
   * @param window_length length of the rolling hash
   * @return size of hash table
   */
  unsigned Initialize(unsigned size_of_table, uint32 window_length);
  /**
   * Initializes the value of rolling hash with the hash value of
   * [*from, *(from +len) ).
   *
   * @param from pointer to the start of sequence
   * @param len length of sequence
   * @return hash value of sequence [*from, *(from +len) )
   */
  int64 InitValue(const byte *from, unsigned len);
  /**
   * Get the expected size of the hash table.
   *
   * @return expected size of the hash table
   */
  unsigned Size() {
    return mask_+1;
  }
  /**
   * Updates the hash value from the hash of [old_val, new_val) to
   * (old_val, new_val].
   *
   * @return updated hash value
   */
  int64 Update(byte old_val, byte new_val);

 private:
  /**<Hash value of previous sequence*/
  int64 prev_hash_;
  /**<Parameter used in calculation of hash function*/
  int64 c_; 
  /**<Mask used in the calculation of hash function. Hash values will
     be in range [0,mask_] */
  int64 mask_; 

};

} // namespace hash_functions

} //namespace long_sequences
} //namespace bwtc

#endif
