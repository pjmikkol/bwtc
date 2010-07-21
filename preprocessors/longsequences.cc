#include <cassert>
#include <cstdlib>

#include <vector>

#include "longsequences.h"
#include "../globaldefs.h"

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

/* Finds prime which is greater or equal to given n */
int FindPrimeGeq(int n) {
  if (n%2 == 0) n += 1;
  while (!IsPrime(n)) n += 2;
  return n;
}

/* Hashtable based structure for holding the frequency information about the *
 * substrings encountered. Calculating the hash values is based on Karp-Rabin*
 * pattern-matching algorithm.                                               *
 *                                                                           *
 * */
class SequenceTable {

  struct seq {
    seq() : count(0), overflow(0), position(0) {}
    seq(unsigned c, int o, uint64 p) : count(c), overflow(o), position(p) {}
    
    unsigned count; /* How many times we have encountered the string */
    int overflow; /* Index of overflow sequence in over_ -vector */
    uint64 position; /* Position in actual text */
  };

 public:
  explicit SequenceTable(unsigned size, unsigned block_length, byte *data)
      : data_(data), block_length_(block_length), prev_hash_(0),
        prev_position_(0)
  {
    size_ = FindPrimeGeq(size);
    table_ = new seq[size_];
    c_ = 1;
    for(unsigned i = 0; i < block_length_ - 1; ++i) {
      c_ <<= 8;
      c_ %= size_;
    }
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
    table_[prev_hash_] = seq(1, -1, 0);
  }

  int64 Search(uint64 pos) {
    UpdateHashValue(pos);
    seq *candidate = &table_[prev_hash_];
    while(candidate->count > 0 && candidate->overflow >= 0) {
      if (StrEq(pos, candidate->position)) {
        return candidate->position;
      }
      else {
        candidate = &over_[candidate->overflow];
      }
    }
    return pos;    
  }

  int64 Insert(uint64 pos) {
    UpdateHashValue(pos);
    seq *candidate = &table_[prev_hash_];
    while(candidate->count > 0 && candidate->overflow >= 0) {
      if (StrEq(pos, candidate->position)) {
        ++candidate->count;
        return candidate->position;
      }
      else {
        candidate = &over_[candidate->overflow];
      }
    }
    candidate->overflow = over_.size();
    over_.push_back(seq(1, -1, pos));
    return pos;
  }

 private:
  bool StrEq(uint64 i, uint64 j) {
    for(unsigned k = 0; k < block_length_; ++k) {
      if(data_[i++] != data_[j++]) return false;
    }
    return true;
  }
  
  void UpdateHashValue(uint64 pos) {
    assert(pos > 0);
    prev_hash_ -= c_*data_[pos - 1];
    prev_hash_ <<= 8;
    prev_hash_ += data_[pos + block_length_ - 1];
    prev_hash_ %= size_;
    if (prev_hash_ < 0) prev_hash_ += size_;
    assert(prev_hash_ >= 0);
  }
  
  byte *data_; /* Source of strings stored to table */
  unsigned block_length_; /* Length of stored strings */
  int64 prev_hash_; /* Previous value of the hash function */
  uint64 prev_position_; /* Position of the previous value in data_ */
  int64 c_; /* parameter used in updating the hash values */ 
  unsigned size_; /* Size of actual hashtable */
  seq *table_; /* hashtable */
  std::vector<seq> over_; /* overflow- lists*/
};

} //namespace long_sequences

/* Uses modification of Bentley-McIlroy algorithm which is based on Karp-Rabin *
 * pattern-matching algorithm.  */
uint64 CompressSequences(byte *from, uint64 length, int memory_constraint)
{
  using namespace long_sequences;

  assert(length > 0);
  assert(from);
  /****************************************************************************
   * We have to store at most length/block_size values to SequenceTable. Each *
   * entry needs 16 bytes of memory.                                          *
   * Let 'a' = (length/block_size)/size_of_table and let c be the coefficient *
   * for memory usage (we can use c*length bytes). Then we can be sure that   *
   *               'a' >= 16/(c*m)                                            *
   * It also holds that given the memory requirements                         *
   *               'a' <= 16/(c*m - 16)                                       *
   ****************************************************************************/
  /* m = 16, c = memory_constraint */
  unsigned block_length = 16;
  assert(length > block_length);
  assert(memory_constraint > 1);
  /* At the moment ensures that the max memory SequenceTable uses is
   * memory_constraint*length bytes */
  SequenceTable seq_table((length/(16*block_length))*
                          (memory_constraint*block_length - 16),
                          block_length, from);
  seq_table.Initialize();
  for(uint64 i = 1; i <= length - block_length; ++i) {
    uint64 prev_occ;
    if( i % block_length == 0) {
      prev_occ = seq_table.Insert(i);
      if(prev_occ != i) {
        /* Check the long sequence */
      }
    }
    else {
      prev_occ = seq_table.Search(i);
      if(prev_occ != i) {
        /* Check the long sequence */        
      }
    }
  }

  return length;
}


} // namespace bwtc
