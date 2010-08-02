#include <cassert>
#include <cstdlib>

#include <list>
#include <utility>
#include <vector>

#include "longsequences.h"
#include "preprocessor.h"
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
  assert(n > 3);
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
      : data_(data), block_length_(block_length), prev_hash_(0)
  {
    //TODO: what if given size too small
    size_ = FindPrimeLeq((memory_constraint*size_of_data)/16 -
                         size_of_data/(2*block_length_));
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

  /* Initializes the structure */
  void Initialize() {
    prev_hash_ = 0;
    for(unsigned i = 0; i < block_length_; ++i) {
      prev_hash_ <<= 8;
      prev_hash_ += data_[i];
      prev_hash_ %= size_;
    }
    table_[prev_hash_] = seq(1, prev_hash_, 0);
  }

  void JumpToPos(uint64 pos) {
    prev_hash_ = 0;
    for(uint64 i = pos; i < pos + block_length_; ++i) {
      prev_hash_ <<= 8;
      prev_hash_ += data_[i];
      prev_hash_ %= size_;
    }
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

  void DebugPrint() {
    for(unsigned i = 0; i < size_; ++i) {
      if(table_[i].count > 0 && (table_[i].overflow & 0x80000000) == 0 ) {
        unsigned j = i;
        do {
          std::cout << i << "\t" << table_[i].position << "\t"
                    <<  table_[i].count << "\n";
          i = table_[i].overflow & 0x7fffffff;
        } while (i != j);
        std::cout << "\n";
      }
    }
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
  int64 c_; /* parameter used in updating the hash values */ 
  unsigned size_; /* size of the whole array */
  unsigned h_size_; /* size of the actual hashtable */
  seq *table_; 
  unsigned next_overflow_; /* next free spot to store overflow-entry */
};

template <typename T>
class CircularBuffer {
 public:
  CircularBuffer(unsigned size) : buffer_(size), head_(size-1), size_(size) {}

  void Insert(T elem) {
    buffer_[head_++] = elem;
    head_ %= size_;
  }

  T Head() const {
    return buffer_[head_];
  }

  const T& operator[](int i) const {
    return buffer_[(head_ + i) % size_];
  }

  void Fill(T value) {
    std::fill(buffer_.begin(), buffer_.end(), value);
  }

  void SetHeadAndRest(T head, T rest) {
    std::fill(buffer_.begin(), buffer_.end(), rest);
    buffer_[head_] = head;
  }

 private:
  std::vector<T> buffer_;
  unsigned head_;
  unsigned size_;
};

/* Data structure for holding long sequences. At the moment searching is *
 * done in simple brute-force style. */
class LongSequences {
  struct long_seq {
    long_seq() : position(0), length(0), count(0) {}
    long_seq(unsigned c, int l, uint64 p) :
        position(p), length(l), count(c) {}

    bool operator<(const long_seq& ll) const {
      return position < ll.position;
    }
    uint64 position;
    int length;
    unsigned count;
  };

 public:
  explicit LongSequences(byte *data, int threshold)
      : data_(data), threshold_(threshold) {}
  
  int64 Update(uint64 pos, int length, unsigned amount) {
    assert(length >= threshold_);
    for(std::list<long_seq>::iterator it = sequences_.begin();
        it != sequences_.end(); ++it)
    {
      std::pair<uint64, int> match = Match(pos, length, *it);
      if (match.second >= threshold_) {
        /* In some cases we split the match */
        int l_leftover = match.first - pos;
        int r_leftover = length - match.second;
        if (l_leftover > 0) r_leftover -= l_leftover;
        if( l_leftover < threshold_ && r_leftover < threshold_) {
          int left_len = match.first - it->position;
          int right_len = it->length - left_len - match.second;
          if(left_len >= threshold_ || right_len >= threshold_) {
            long_seq old = *it;
            sequences_.erase(it);
            if(left_len >= threshold_)
              Update(old.position, left_len, old.count);
            Update(match.first, match.second, old.count + amount);
            if(right_len >= threshold_)
              return Update(match.first + match.second, right_len, old.count);
          } else {
            it->count += amount;
            it->position = match.first;
            it->length = match.second;
          }
        } else {
          it->count += amount;
          it->position = match.first;
          it->length = match.second;
          if(l_leftover >= threshold_) Update(pos, l_leftover, amount);
          assert(r_leftover < length);
          if(r_leftover >= threshold_)
            return Update(match.second + match.first, r_leftover, amount);
        }
        return match.first + match.second;
      }
    }
    if(amount == 1) ++amount;
    sequences_.push_back(long_seq(amount, length, pos));
    return pos + length;
  }

  void DebugPrint() {
    unsigned s = 0;
    for(std::list<long_seq>::iterator it = sequences_.begin();
        it != sequences_.end(); ++it)
    {
      std::cout << "Frequency: " << it->count << " Length: "
                << it->length <<"\n";
      for(uint64 i = it->position; i < it->position + it->length; ++i)
        std::cout << data_[i];
      std::cout << "\n-----------------------\n";
      ++s;
    }
    std::cout << "Total size of table: " << s << "\n";
  }
    
 private:
  /* Returns pair <length,starting position> for the match. Uses brute-force  *
   * strategy. Returns first match which has length is greater than threshold_*/
  std::pair<uint64, int> Match(uint64 pos, int length,
                                    const long_seq& sequence)
  {
    if(length < sequence.length) {
      for(int i = 0; i < sequence.length - threshold_; ++i) {
        for(int j = 0; j < length - threshold_; ++j) {
          int k = i, l = j;
          while(l < length && data_[pos + l] == data_[sequence.position + k]) {
            ++l; ++k;
          }
          if(l - j >= threshold_)
            return std::pair<uint64, int>(sequence.position + i, l - j);
        }
      }
    } else {
      for(int i = 0; i < length - threshold_; ++i) {
        for(int j = 0; j < sequence.length - threshold_; ++j) {
          int k = i, l = j;
          while(l < sequence.length &&
                data_[pos + k] == data_[sequence.position + l])
          {
            ++l; ++k;
          }
          if( l - j >= threshold_)
            return std::pair<uint64, int>(sequence.position + j, l - j);
        }
      }
    }
    return std::pair<uint64, int>(0,0);
  }

  byte *data_;
  std::list<long_seq> sequences_;
  int threshold_; /* Threshold length for sequence to be long enough */
};

/* Calculates frequencies from range [source_begin, source_end) */
void CalculateFrequencies(uint64 *target, byte *source_begin, byte *source_end) {
  while(source_begin != source_end) ++target[*source_begin++];
}

void DetectSequences(byte *from, uint64 length, int memory_constraint,
                     unsigned block_length, int threshold, uint64 *freqs,
                     LongSequences *long_seqs)
{
  /* Initialize the data structures */
  SequenceTable seq_table(length, block_length, from, memory_constraint);
  seq_table.Initialize();
  CircularBuffer<bool> buffer(block_length);
  buffer.Fill(false);
  seq_table.JumpToPos(block_length - 1);
  CalculateFrequencies(freqs, from, from + block_length);
  
  int matches = 0, longs = 0;
  int64 long_end = 0;
  for(uint64 i = block_length; i <= length - block_length; ++i) {
    ++freqs[from[i]];
    uint64 prev_occ;
    if( i % block_length == 0) prev_occ = seq_table.Insert(i);
    else prev_occ = seq_table.Search(i);
    if(prev_occ != i) {
      ++matches;
      if( buffer.Head() ) {
        /* Check for the long sequence */
        uint64 k = prev_occ;
        int64 start_pos = i;
        while(k > 0 && start_pos >= long_end /*start_pos > k + block_length */
              && from[k] == from[start_pos]) {
          --k; --start_pos;
        }
        if( start_pos < 0 || from[k] != from[start_pos])
            /*(k == 0 && from[k] != from[start_pos]) ||
              (k != 0 && start_pos != k + block_length))*/
          ++start_pos;
        uint64 end_pos = i + block_length;
        k = prev_occ + block_length;
        while(end_pos < length && k < start_pos && from[k] == from[end_pos]) {
          ++k; ++end_pos;
        }
        if((end_pos == length && from[k] != from[end_pos]) ||
           (end_pos != length && k < start_pos))
          --end_pos;
        int seq_len = end_pos - start_pos + 1;
        if(seq_len >= threshold) {
          ++longs;
          long_end = long_seqs->Update(start_pos, seq_len, 1);
          seq_table.JumpToPos(end_pos);
          buffer.Fill(false);
          CalculateFrequencies(freqs, from + i + 1, from + end_pos + 1);
          i = end_pos;
          continue;
        }
      }
      /* Which one of the alternative ways is better? 
      buffer.SetHeadAndRest(true, false);
      CalculateFrequencies(freqs, from + i + 1, from + i + block_length);
      i += block_length - 1;
      seq_table.JumpToPos(i);*/
      buffer.Insert(true);
    } else
      buffer.Insert(false);
  }
  //seq_table.DebugPrint();
  long_seqs->DebugPrint();
  std::cout << "Total matches: " << matches << "\n";
  std::cout << "Long matches: " << longs << "\n";
  
}

} //namespace long_sequences

/* Uses modification of Bentley-McIlroy algorithm which is based on Karp-Rabin *
 * pattern-matching algorithm.  */
uint64 CompressSequences(byte *from, uint64 length, int memory_constraint,
                         unsigned block_length, int threshold)
{
  using namespace long_sequences;
  assert(length > 0);
  assert(from);
  assert(length > block_length);
  assert(memory_constraint > 1);
  LongSequences long_seqs(from, threshold);
  uint64 frequencies[256] = {0};

  //TODO: figure how the function will return information about sequences
  //      of length block_length
  DetectSequences(from, length, memory_constraint, block_length, threshold,
                  frequencies, &long_seqs);


  return length;
}


} // namespace bwtc
