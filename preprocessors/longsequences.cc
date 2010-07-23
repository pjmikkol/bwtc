#include <cassert>
#include <cstdlib>

#include <vector>

#include "longsequences.h"
#include "../globaldefs.h"

#include <iostream>

namespace bwtc {

namespace long_sequences {

namespace {
/* Utilities for finding prime numbers. Test we use here is deterministic *
 * variant of Miller-Rabin primality test. */

int mulm(int a, int b, int m)
{
  return static_cast<int64>(a)*b % m;
}

template<class T>
T powm(T b, T e, T m)
{
  if (e == 0) return 1;
  T x = powm(b, e>>1, m);
  T xx = mulm(x, x, m);
  return (e & 1) ? mulm(b, xx, m) : xx;
}

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


} // empty namespace

bool IsPrime(int n)
{
  assert(n % 2 == 1);
  if ( n < 1373653) return (!Witness(2,n) && !Witness(3,n));
  else if ( n < 9080191 ) return (!Witness(31,n) && !Witness(73,n));
  /* n < 4,759,123,141 */
  else return (!Witness(2,n) && !Witness(7,n) && !Witness(61,n));
}

/* Finds prime which is lesser or equal to given n */
int FindPrimeLeq(int n) {
  assert(n > 4);
  if (n%2 == 0) n -= 1;
  while (!IsPrime(n)) n -= 2;
  return n;
}

/* Hashtable based structure for holding the frequency information about the *
 * substrings encountered. Calculating the hash values is based on Karp-Rabin*
 * pattern-matching algorithm.                                               *
 *                                                                           *
 * Same object is used also for finding the most frequent strings.           */
class SequenceTable {

  struct seq {
    seq() : count(0), overflow(0), position(0) {}
    seq(unsigned c, unsigned o, uint64 p) :
        count(c), overflow(o), position(p) {}
    
    unsigned count; /* How many times we have encountered the string */
    unsigned overflow; /* Index of overflow sequence in over_ -vector */
    uint64 position; /* Position in actual text */
  };

 public:
  /****************************************************************************
   * Let D be the number of bytes in the data and let B be the size of block. * 
   * We have to store at most D/B values to SequenceTable. Single entry takes *
   * 16 bytes. If we allocate 0.5*D/B places for overflow-part of the array   *
   * then to satisfy memory constraint we can set the actual size of the      *
   * hashtable to (c*D)/16 - D/(2*m) where c is how many bytes we can use for *
   * single byte in input.                                                    *
   ****************************************************************************/
  explicit SequenceTable(unsigned size_of_data, unsigned block_length,
                         byte *data, unsigned memory_constraint)
      : data_(data), block_length_(block_length), prev_hash_(0),
        prev_position_(0)
  {
    //TODO: what if given size too small
    size_ = FindPrimeLeq((memory_constraint*size_of_data)/16.0 -
                         size_of_data/(2.0*block_length_));
    h_size_ = size_ + size_of_data/(2*block_length_);
    table_ = new seq[h_size_];
    next_overflow_ = h_size_ - 1;
    c_ = 1;
    for(unsigned i = 0; i < block_length_ - 1; ++i) {
      c_ <<= 8;
      c_ %= size_;
    }
    assert((size_ & 0x80000000) == 0);
    assert((h_size_ & 0x80000000) == 0);
  }

  ~SequenceTable() {
    delete [] table_;
  }

  /* Initializes the structure and read*/
  void Initialize() {
    prev_hash_ = 0;
    for(unsigned i = 0; i < block_length_; ++i) {
      prev_hash_ <<= 8;
      prev_hash_ += data_[i];
      prev_hash_ %= size_;
    }
    prev_position_ = 0;
    table_[prev_hash_] = seq(1, prev_hash_, 0);
  }

  int64 Search(uint64 pos) {
    UpdateHashValue(pos);
    seq *head = &table_[prev_hash_];
    if(head->count == 0) return pos;
    seq *curr = head;
    do {
      if (StrEq(pos, curr->position)) {
        ++curr->count;
        return curr->position;
      }
      curr = &table_[curr->overflow & 0x7fffffff];
    } while(curr != head);
    return pos;
  }

  int64 Insert(uint64 pos) {
    UpdateHashValue(pos);
    seq *head = &table_[prev_hash_];
    if (head->count == 0) {
      table_[prev_hash_] = seq(1, prev_hash_, pos);
      return pos;
    }
    else if ((head->overflow & 0x80000000) == 0) {
      /* There are already entries with the same hash value */
      seq *curr = head;
      do {
        if (StrEq(pos, curr->position)) {
          ++curr->count;
          return curr->position;
        }
        curr = &table_[curr->overflow & 0x7fffffff];
      } while(curr != head);
      
      table_[next_overflow_] = seq(1, 0x80000000 | head->overflow, pos);
      head->overflow =  next_overflow_;
    } else {
      seq *curr = head;
      while ((curr->overflow & 0x7fffffff) != prev_hash_) {
        curr = &table_[curr->overflow & 0x7fffffff];
      }
      table_[next_overflow_] = *head;
      curr->overflow = (curr->overflow & 0x80000000) | next_overflow_;
      table_[prev_hash_] = seq(1, prev_hash_, pos);
    }
    FindNextOverflowPos();
    return pos;
  }

 private:
  void FindNextOverflowPos() {
    while(table_[next_overflow_].count > 0 ) --next_overflow_;
  }
  
  bool StrEq(uint64 i, uint64 j) {
    for(unsigned k = 0; k < block_length_; ++k) {
      if(data_[i++] != data_[j++]) return false;
    }
    return true;
  }
  
  void UpdateHashValue(uint64 pos) {
    assert(pos > 0);
    prev_hash_ -= c_*data_[pos - 1];
    prev_hash_ <<= 8; /* Multiply with the size of alphabet */
    prev_hash_ += data_[pos + block_length_ - 1];
    prev_hash_ %= size_;
    if (prev_hash_ < 0) prev_hash_ += size_;
    assert(prev_hash_ >= 0);
  }
  /* We use one array for holding the hashtable and its overflow-lists. *
   * Range [0,h_size_) is reserved for hashtable. Initially range       *
   * [h_size_, size_) is for overflow-lists. If there are more than     *
   * size_ - h_size_ entries in overflow-lists we start to use actual   *
   * hashtable for overflow-entries. First bit of entry's overflow-field*
   * tells whether or not this entry belongs to overflow-list           */
  byte *data_; /* Source of strings stored to table */
  unsigned block_length_; /* Length of stored strings */
  int64 prev_hash_; /* Previous value of the hash function */
  uint64 prev_position_; /* Position of the previous value in data_ */
  int64 c_; /* parameter used in updating the hash values */ 
  unsigned size_; /* size of the whole array */
  unsigned h_size_; /* size of the actual hashtable */
  seq *table_; 
  unsigned next_overflow_; /* next free spot to store overflow-entry */
};

} //namespace long_sequences

/* Uses modification of Bentley-McIlroy algorithm which is based on Karp-Rabin *
 * pattern-matching algorithm.  */
uint64 CompressSequences(byte *from, uint64 length, int memory_constraint)
{
  using namespace long_sequences;

  assert(length > 0);
  assert(from);
  unsigned block_length = 16;
  assert(length > block_length);
  assert(memory_constraint > 1);
  SequenceTable seq_table(length, block_length, from, memory_constraint);
  seq_table.Initialize();
  for(uint64 i = 1; i <= length - block_length; ++i) {
    uint64 prev_occ;
    if( i % block_length == 0) {
      prev_occ = seq_table.Insert(i);
      if(prev_occ != i) {
        /* Check for the long sequence */
      }
    }
    else {
      prev_occ = seq_table.Search(i);
      if(prev_occ != i) {
        /* Check for the long sequence */        
      }
    }
  }
  




  return length;
}


} // namespace bwtc
