/**************************************************************************
 *  Copyright 2010, Pekka Mikkola, pjmikkol (at) cs.helsinki.fi           *
 *                                                                        *
 *  This file is part of bwtc.                                            *
 *                                                                        *
 *  bwtc is free software: you can redistribute it and/or modify          *
 *  it under the terms of the GNU General Public License as published by  *
 *  the Free Software Foundation, either version 3 of the License, or     *
 *  (at your option) any later version.                                   *
 *                                                                        *
 *  bwtc is distributed in the hope that it will be useful,               *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *  GNU General Public License for more details.                          *
 *                                                                        *
 *  You should have received a copy of the GNU General Public License     *
 *  along with bwtc.  If not, see <http://www.gnu.org/licenses/>.         *
 **************************************************************************/

#include <cassert>
#include <cstdlib>

#include <algorithm>
#include <list>
#include <utility>
#include <vector>

#include "longsequences.h"
#include "preprocessor.h"
#include "../coders.h"
#include "../globaldefs.h"
#include "../utils.h"

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
  if ((n&1) == 0) n -= 1;
  while (!IsPrime(n)) n -= 2;
  return n;
}

/* Base class for computing different rolling hashes */
class Hasher {
 public:
  virtual ~Hasher() {};
  /* Initialize function computes the size of actual hash_table and space *
   * reserved for overflow-lists. Their sum is returned.                  */
  virtual unsigned Initialize(unsigned max_values, int64 bytes_to_use,
                              unsigned size_of_elem, unsigned w_length) = 0;
  virtual int64 InitValue(const byte *from, unsigned length) = 0;
  virtual int64 Update(byte old_val, byte new_val) = 0;

};
  
class PrimeHasher : public Hasher {
 public:
  PrimeHasher() : prev_hash_(0) {}
  virtual ~PrimeHasher() {}
  virtual unsigned Initialize(unsigned max_values, int64 bytes_to_use,
                              unsigned size_of_elem, unsigned window_length)
  {
    q_ = FindPrimeLeq(bytes_to_use/size_of_elem - max_values/2);
    c_ = 1;
    for(unsigned i = 0; i < window_length - 1; ++i) {
      c_ <<= 8;
      c_ %= q_;
    }
    assert(q_ <= bytes_to_use/size_of_elem);
    //return bytes_to_use/size_of_elem;
    return q_ + max_values/2;
  }

  virtual int64 InitValue(const byte *from, unsigned len) {
    prev_hash_ = 0;
    for(unsigned i = 0; i < len; ++i)
      prev_hash_ = ((prev_hash_ << 8) + *from++) % q_;
    return prev_hash_;
  }

  virtual int64 Update(byte old_val, byte new_val) {
    prev_hash_ = (((prev_hash_ - old_val*c_) << 8) + new_val) % q_;
    if (prev_hash_ < 0) prev_hash_ = -prev_hash_;
    return prev_hash_;
  }

 private:
  int64 prev_hash_;
  int64 c_; /* parameter used in calclation of hash values */
  int64 q_; /* prime, hash values will be in range [0..q_) */
};

class FastHasher : public Hasher {
 public:
  FastHasher() : prev_hash_(0) {}
  virtual ~FastHasher() {}
  virtual unsigned Initialize(unsigned max_values, int64 bytes_to_use,
                              unsigned size_of_elem, unsigned window_length)
  {
    //TODO: Put more care for choosing the mask_ and size of hash table
    mask_ = MostSignificantBit(max_values) - 1;
    c_ = 1;
    for(unsigned i = 0; i < window_length - 1; ++i)
      c_ = (c_*257) & mask_;
    assert(mask_ <= bytes_to_use/size_of_elem);
    return bytes_to_use/size_of_elem;
  }

  virtual int64 InitValue(const byte *from, unsigned len) {
    prev_hash_ = 0;
    for(unsigned i = 0; i < len; ++i)
      prev_hash_ = ((prev_hash_*257) + *from++) & mask_;
    return prev_hash_;
  }

  virtual int64 Update(byte old_val, byte new_val) {
    prev_hash_ = (((prev_hash_ - old_val*c_)*257) + new_val) & mask_;
    return prev_hash_;
  }

 private:
  int64 prev_hash_;
  int64 c_; /* parameter used in calclation of hash values */
  int64 mask_; /* prime */
};


/******************************************************************************
 * Hashtable based structure for holding the frequency information about the  *
 * substrings encountered. Calculating the hash values is based on Karp-Rabin *
 * pattern-matching algorithm.                                                *
 *                                                                            *
 * Same object is used also for finding the most frequent strings.            *
 ******************************************************************************/
class SequenceTable {

  struct seq {
    seq() : count(0), overflow(0), position(0) {}
    seq(unsigned c, unsigned o, uint64 p) :
        count(c), overflow(o), position(p) {}
    
    unsigned count; /* How many times we have encountered the string */
    unsigned overflow; /* Index of overflow sequence in over_ -vector */
    unsigned position; /* Position in actual text */
  };

 public:
  /****************************************************************************
   * Let D be the number of bytes in the data and let B be the size of block. * 
   * We have to store at most D/B values to SequenceTable. Single entry takes *
   * 16 bytes. If we allocate 0.5*D/B places for overflow-part of the array   *
   * then to satisfy memory constraint we can set the actual size of the      *
   * hashtable to (c*D)/s - D/(2*m) where c is how many bytes we can use for  *
   * single byte in input and s is the size of single seq-struct.             *
   *                                                                          *
   * Parameter hasher decides which hash function to use.                     *
   ****************************************************************************/
  explicit SequenceTable(uint64 size_of_data, unsigned block_length,
                         byte *data, unsigned memory_constraint,
                         char hasher = 'p')
      : data_(data), block_length_(block_length), prev_hash_(0), h_(0)
  {
    //TODO: what if given size too small
    switch(hasher) {
      case 'p':
        h_ = new PrimeHasher();
        break;
      default:
        h_ = new FastHasher();
        break;
    }
    h_size_ = h_->Initialize(size_of_data/block_length_,
                             memory_constraint*size_of_data, sizeof(seq),
                             block_length);
    assert(h_size_ >= size_of_data/block_length);
    assert(h_size_*sizeof(seq) > size_of_data/block_length);
    table_ = new seq[h_size_];
    next_overflow_ = h_size_ - 1;
    assert((h_size_ & 0x80000000) == 0);
  }

  ~SequenceTable() {
    delete [] table_;
  }

  /* Initializes the structure */
  void Initialize() {
    prev_hash_ = h_->InitValue(data_, block_length_);
    table_[prev_hash_] = seq(1, prev_hash_, 0);
  }

  void JumpToPos(uint64 pos) {
    prev_hash_ = h_->InitValue(data_ + pos, block_length_);
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

  void Periods(std::vector<long_seq> *periods) {
    /*
    for(unsigned i = 0; i < h_size_; ++i) {
      if(table_[i].count > 0) {
        periods->push_back(
            long_seq(table_[i].count, block_length_, table_[i].position));
        }
    }
    */
    Border<byte> b(block_length_);
    for(unsigned i = 0; i < h_size_; ++i) {
      if(table_[i].count > 0) {
        unsigned border = b(data_ + table_[i].position);
        if (border >= block_length_ - 1) {
          periods->push_back(
              long_seq(table_[i].count, 2, table_[i].position));
        } else {
          periods->push_back(
              long_seq(table_[i].count, block_length_ - border,
                       table_[i].position));
        }
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
    prev_hash_ = h_->Update(data_[pos - 1], data_[pos + block_length_ - 1]);
  }
  /***********************************************************************
   * We use one array for holding the hashtable and its overflow-lists.  *
   * Range [0,x) is reserved for hashtable. Initially range              *
   * [x, h_size_) is for overflow-lists. If there are more than          *
   * h_size_ - x entries in overflow-lists we start to use actual        *
   * hashtable for overflow-entries. Only Hasher knows value of x.       *
   * First bit of entry's overflow-field tells whether or not this entry *
   * belongs to overflow-list                                            *
   ***********************************************************************/
  byte *data_; /* Source of strings stored to table */
  unsigned block_length_; /* Length of stored strings */
  int64 prev_hash_; /* Previous value of the hash function */
  int64 c_; /* parameter used in updating the hash values */ 
  unsigned h_size_; /* size of the whole array */
  seq *table_; 
  unsigned next_overflow_; /* next free spot to store overflow-entry */
  Hasher *h_;
};

template <typename T>
class CircularBuffer {
 public:
  CircularBuffer(unsigned size) : head_(size-1), size_(size) {
    buffer_ = new T[size_];
  }

  ~CircularBuffer() {
    delete [] buffer_;
  }

  void Insert(T elem) {
    buffer_[head_++] = elem;
    if(size_ == head_) head_ = 0;
  }

  T Head() const {
    return buffer_[head_];
  }

  const T& operator[](int i) const {
    return buffer_[(head_ + i) % size_];
  }

  void Fill(T value) {
    std::fill(buffer_, buffer_ + size_, value);
  }

  void SetHeadAndRest(T head, T rest) {
    std::fill(buffer_, buffer_ + size_, rest);
    buffer_[head_] = head;
  }

 private:
  unsigned head_;
  unsigned size_;
  T *buffer_;
};

/* Replacement struct is used in lists which are inspected when writing *
 * replcements of sequences. */
struct replacement {
  explicit replacement(int64 loc, unsigned len, byte s) :
      location(loc), length(len), rank(s) {}

  int64 location;
  unsigned length;
  byte rank;
};

/* We do the comparison of the srings backwards, because it seems to be *
 * faster then forwards. */
bool SeqEq(const long_seq& s1, const long_seq& s2, byte *data) {
  if (s1.length != s2.length) return false;
  byte *f = data + s1.position + s1.length - 1;
  byte *g = data + s2.position + s1.length - 1;
  const byte *start = data + s1.position - 1;
  while(f != start) if(*f-- != *g--) return false;
  /*  for(int i = 0; i < s1.length; ++i)
    if(data[s1.position + i] != data[s2.position + i])
      return false;
      }*/
  return true;
}

bool SeqEq(const replacement& s, uint64 position, byte *data,
           uint64 len_of_data)
{
  if (s.length + position > len_of_data) return false;
  /*
  if (s.length + position > len_of_data) return false;
  for(unsigned i = 0; i < s.length; ++i) {
    if(data[s.location + i] != data[position + i])
      return false;
      }*/
  //assert(data[s.location] == data[position]);
  //assert(data[s.location + 1] == data[position + 1]);
  byte *start = data + s.location + 1;
  byte *f = start + s.length - 2;
  byte *g = data + position + s.length - 1;
  while(f != start) if(*f-- != *g--) return false;
  return true;
}

/* Data structure for holding long sequences. At the moment searching is *
 * done in simple brute-force style. */
class LongSequences {
  struct match {
    match() : match_pos(0), orig_pos(0), length(0) {}
    match(uint64 m_pos, uint64 o_pos, unsigned l) :
        match_pos(m_pos), orig_pos(o_pos), length(l) {}
    uint64 match_pos;
    uint64 orig_pos;
    unsigned length;
  };

  struct cmp_long_seq {
    cmp_long_seq(byte *s) : source(s) {}

    bool operator()(const long_seq& s1, const long_seq& s2) const {
      int min_len = std::min(s1.length, s2.length);
      for(int j = 0; j < min_len; ++j) {
        if(source[s1.position + j] != source[s2.position + j]) {
          return source[s1.position + j] < source[s2.position + j];
        }
      }
      return s1.length < s2.length;
    }

    byte *source;
  };

 public:
  explicit LongSequences(byte *data, int threshold)
      : data_(data), threshold_(threshold) {}
  
  /* Inefficient implementation, byt should be only callled once during the *
   * algorithm. */
  void MergeDuplicatesAndSort() {
    size_t s = sequences_.size();
    if(s <= 1) return;
    std::vector<long_seq> temp;
    temp.resize(s);
    std::copy(sequences_.begin(), sequences_.end(), temp.begin());
    std::sort(temp.begin(), temp.end(), cmp_long_seq(data_));
    std::vector<long_seq> temp2;
    std::vector<long_seq>::iterator curr = temp.begin();
    std::vector<long_seq>::iterator prev = curr++;
    for(; curr != temp.end(); prev = curr++)
    {
      if(SeqEq(*prev, *curr, data_)) {
        curr->count += prev->count - 1;
      } else  {
        temp2.push_back(*prev);
      }
    }
    temp2.push_back(*prev);
    std::sort(temp2.begin(), temp2.end());
    std::list<long_seq>::iterator lit = sequences_.begin();
    for(curr = temp2.begin(); curr != temp2.end(); ++curr, ++lit) {
      *lit = *curr;
    }
    sequences_.erase(lit, sequences_.end());
  }

  bool Empty() {
    return sequences_.empty();
  }
  
  long_seq Pop() {
    long_seq val = sequences_.back();
    sequences_.pop_back();
    return val;
  }

  int64 Update(uint64 pos, int length, unsigned amount) {
    assert(length >= threshold_);
    for(std::list<long_seq>::iterator it = sequences_.begin();
        it != sequences_.end(); ++it)
    {
      match m = Match(pos, length, *it);
      if (m.length >= static_cast<unsigned>(threshold_) ) {
        /* In some cases we split the match */
        int l_leftover = m.match_pos - pos;
        int r_leftover = length - m.length;
        if(l_leftover > 0) r_leftover -= l_leftover;
        if(l_leftover < threshold_ && r_leftover < threshold_) {
          int left_len = m.orig_pos - it->position;
          int right_len = it->length - left_len - m.length;
          if(left_len >= threshold_ || right_len >= threshold_) {
            long_seq old = *it;
            sequences_.erase(it);
            if(left_len >= threshold_)
              Update(old.position, left_len, old.count);
            if(right_len >= threshold_)
              Update(m.orig_pos + m.length, right_len, old.count);
            return Update(m.match_pos, m.length, old.count + amount);
          } else {
            it->count += amount;
            it->position = m.match_pos;
            it->length = m.length;
          }
        } else {
          it->count += amount;
          it->position = m.match_pos;
          it->length = m.length;
          if(l_leftover >= threshold_) Update(pos, l_leftover, amount);
          assert(r_leftover < length);
          if(r_leftover >= threshold_)
            return Update(m.length + m.match_pos, r_leftover, amount);
        }
        return m.match_pos + m.length;
      }
    }
    if(amount == 1) ++amount;
    sequences_.push_front(long_seq(amount, length, pos));
    return pos + length;
  }

  void DebugPrint() {
    unsigned s = 0;
    for(std::list<long_seq>::iterator it = sequences_.begin();
        it != sequences_.end(); ++it)
    {
      std::clog << "Frequency: " << it->count << " Length: "
                << it->length << " Position: " << it->position << "\n";
      for(uint64 i = it->position; i < it->position + it->length; ++i)
        std::clog << data_[i];
      std::clog << "\n-----------------------\n";
      ++s;
    }
    std::clog << "Total size of table: " << s;
    std::clog << "\n#########################\n";
  }
    
 private:
  /* Returns match-struct where match_pos is the position of the current match*
   * and orig_pos is the starting position of the match stored in structure.  *
   * Uses brute-force  strategy and returns first match which length is       *
   * greater than threshold_*/
  match Match(uint64 pos, int length, const long_seq& sequence)
  {
    if(length < sequence.length) {
      for(int i = 0; i < sequence.length - threshold_; ++i) {
        for(int j = 0; j < length - threshold_; ++j) {
          int k = i, l = j;
          while(l < length && data_[pos + l] == data_[sequence.position + k]) {
            ++l; ++k;
          }
          if(l - j >= threshold_)
            return match(pos + j, sequence.position + i, l - j);
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
            return match(pos + i, sequence.position + j, l - j);
        }
      }
    }
    return match(0, 0, 0);
  }

  byte *data_;
  std::list<long_seq> sequences_;
  int threshold_; /* Threshold length for sequence to be long enough */
};

/* Calculates frequencies from range [source_begin, source_end) */
template <typename T>
void CalculateFrequencies(T *target, byte *source_begin, byte *source_end) {
  while(source_begin != source_end) ++target[*source_begin++];
}

void DetectSequences(byte *from, uint64 length, int memory_constraint,
                     unsigned block_length, int threshold, uint64 *freqs,
                     LongSequences *long_seqs, std::vector<long_seq> *periods)
{
  /* Initialize the data structures */
  SequenceTable seq_table(length, block_length, from, memory_constraint);
  seq_table.Initialize();
  CircularBuffer<byte> buffer(block_length);
  buffer.Fill(false);
  seq_table.JumpToPos(block_length - 1);
  CalculateFrequencies(freqs, from, from + block_length);
  
  int matches = 0;
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
        int64 k = prev_occ;
        int64 start_pos = i;
        assert(start_pos > static_cast<int64>(k));
        while(k > 0 && start_pos >= long_end &&  from[k] == from[start_pos] &&
              start_pos > static_cast<int64>(prev_occ))
              /*start_pos > k + block_length */
        {
          --k; --start_pos;
        }
        if(start_pos < long_end || start_pos == static_cast<int64>(prev_occ) ||
           from[k] != from[start_pos])
          ++start_pos;
        uint64 end_pos = i;
        k = prev_occ;
        while(end_pos < length && k < start_pos && from[k] == from[end_pos]) {
          ++k; ++end_pos;
        }
        if(end_pos == length || from[k] != from[end_pos] || k == start_pos)
          --end_pos;
        int seq_len = end_pos - start_pos + 1;
        if(seq_len >= threshold) {
          long_end = long_seqs->Update(start_pos, seq_len, 1);
          seq_table.JumpToPos(end_pos);
          buffer.Fill(false);
          CalculateFrequencies(freqs, from + i + 1, from + end_pos + 1);
          i = end_pos;
          continue;
        }
      }
      /* Which one of the alternative ways is better? */
      buffer.SetHeadAndRest(true, false);
      CalculateFrequencies(freqs, from + i + 1, from + i + block_length);
      i += block_length - 1;
      seq_table.JumpToPos(i);
      //buffer.Insert(true);
    } else
        buffer.Insert(false);
  }
  //seq_table.DebugPrint();
  if(verbosity > 2)
    std::clog << "Total matches in detecting long sequences: "
              << matches << "\n";
  seq_table.Periods(periods);
}

void DecreaseFreqs(FreqTable *freqs, unsigned *vals, unsigned times) {
  for(unsigned i = 0; i < 256; ++i) {
    if(vals[i]) {
      vals[i] *= times;
      freqs->Decrease(i, vals[i]);
    }
  }
}

bool cmp_long_seq_freq(const long_seq& s1, const long_seq& s2) {
    return s1.count*s1.length > s2.count*s2.length;
}

unsigned DecideReplacements(FreqTable *freqs, std::vector<long_seq> *periods,
                            LongSequences *long_seqs,
                            std::list<replacement> *repl, byte *data)
{
  unsigned seq_freqs[256];
  if(verbosity > 6) long_seqs->DebugPrint();
  long_seqs->MergeDuplicatesAndSort();
  if(verbosity > 4) long_seqs->DebugPrint();
  int j = 0;
  while(!long_seqs->Empty() && j < 255) {
    long_seq s = long_seqs->Pop();
    std::fill(seq_freqs, seq_freqs + 256, 0);
    /* We can't update the values because the instances may overlap */
    //CalculateFrequencies(seq_freqs, data + s.position,
    //                   data + s.position + s.length);
    //DecreaseFreqs(freqs, seq_freqs, s.count);
    repl[(data[s.position] << 8) + data[s.position + 1]]
        .push_back(replacement(s.position, s.length, j++));
  }
  if (periods->size() < static_cast<unsigned>(255-j))
    std::sort(periods->begin(), periods->end(), cmp_long_seq_freq);
  else
    std::partial_sort(periods->begin(), periods->begin() + (255-j),
                      periods->end(), cmp_long_seq_freq);
  unsigned i = 0;
  long_seq curr;
  while(j < 255 && i < periods->size()) {
    curr = (*periods)[i];
    if((*freqs)[j] + curr.count >= (curr.count-1)*curr.length) {
      break;
    } else {
      repl[(data[curr.position] << 8) | data[curr.position + 1]].
          push_back(replacement(curr.position, curr.length, j++));
      ++i;
    }
  }
  if(verbosity > 2)
    std::clog << "Replacing " << (j - i) << " sequences longer than "
              << "threshold and " << i << " shorter sequences.\n";

  /* Ensuring that the escape bytes will be written */
  if(j <= 0) return 0xF000;
  if((*freqs)[j-1] > 0) {
    unsigned escape_index = j;
    while(j >= 0 && (*freqs)[j] > 0) {
      unsigned key = freqs->Key(j) << 8;
      for(unsigned k = 0; k < 256; ++k) {
        repl[key | k].push_back(replacement(-1,1,escape_index));
      }
      --j;
    }
    if(verbosity > 2) std::clog << "Using escape byte.\n";
    return escape_index;
  }
  if(verbosity > 2) std::clog << "Not using escape byte.\n";
  return 0xF000;
}

uint64 WriteReplacements(std::list<replacement> *rpls, byte *to, byte *from,
                         uint64 length, byte escape_byte, FreqTable *freqs)
{
  uint64 result_index = 0;
  uint16 pair = static_cast<uint16>(from[0]);
  uint64 i = 1;
  while(1) {
    pair <<= 8;
    pair |= from[i];
    std::list<replacement>::const_iterator it = rpls[pair].begin();
    std::list<replacement>::const_iterator end = rpls[pair].end();
    bool seq_replaced = false;
    while(it != end && it->length > 1) {
      if (SeqEq(*it, i - 1, from, length)) {
        to[result_index++] = freqs->Key(it->rank);
        i += it->length - 1;
        pair = from[i];
        seq_replaced = true;
        break;
      }
      ++it;
    }
    if(it != end && it->length == 1) {
      to[result_index++] = escape_byte;
    }
    if (!seq_replaced) {
      to[result_index++] = from[i-1];
    }
    if( i >= length - 1) {
      assert(i <= length);
      if(i == length) break;
      if(!rpls[from[i] << 8].empty() && rpls[from[i] << 8].back().length == 1)
        to[result_index++] = escape_byte;
      to[result_index++] = from[i];
      break;
    }
    ++i;
  }
  return result_index;
}

} //namespace long_sequences

/* Uses modification of Bentley-McIlroy algorithm which is based on Karp-Rabin *
 * pattern-matching algorithm.  We give priority for replacing the long        *
 * sequences. Sequence is considered long if its length is greate than         *
 * threshold. Handling the long values is quite slow so threshold value should *
 * be quite big. */
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
  std::vector<long_seq> periods;

  DetectSequences(from, length, memory_constraint, block_length, threshold,
                  frequencies, &long_seqs, &periods);
  FreqTable freqs(frequencies);
  std::list<replacement> replacements[65536];
  unsigned escape_index = DecideReplacements(&freqs, &periods, &long_seqs,
                                         replacements, from);
  bool escaping = escape_index < 255;
  byte escape_byte = escaping ? freqs.Key(escape_index) : 0;
  byte *temp = new byte[2*length];
  unsigned position = 0;
  /*************************************************************************
   * Info of the replacements is in a following format:                    *
   * 1)  We write triplets: <s, l, sequence> where sequence is the         *
   *     sequence of length l which is going to be replaced with symbol s. *
   *     Length l is encoded with function utils::PackAndWriteInteger.     *
   * 2a) If we don't do replacements we write only the pair <s, 0>, where  *
   *     s is any byte and 0 is 0-byte.                                    *
   * 2b) Otherwise we signal the end of replacements by writing the symbol *
   *     S which is the same symbol which replaces the previous sequence.  *
   *     (Since we can't replace two sequences with the same symbol, this  *
   *     is ok). After this we write escape symbol if it is in use.        *
   *     Otherwise the S is written again.                                 *
   *************************************************************************/
  unsigned prev_s = 0xF000;
  for(int i = 0; i < 65536; ++i) {
    if(replacements[i].empty()) continue;
    std::list<replacement>::const_iterator it = replacements[i].begin();
    std::list<replacement>::const_iterator end = replacements[i].end();
    while(it != end && it->length > 1) {
      prev_s = freqs.Key(it->rank);
      temp[position++] = prev_s;

      position += utils::PackAndWriteInteger(it->length, temp + position);
      std::copy(&from[it->location], &from[it->location + it->length],
                &temp[position]);
      position += it->length;
      ++it;
    }
  }
  temp[position++] = static_cast<byte>(prev_s & 0xFF);
  if(prev_s == 0xF000) temp[position++] = 0;
  else if(!escaping) temp[position++] = static_cast<byte>(prev_s);
  else temp[position++] = escape_byte;

  uint64 result_length = position;
  result_length += WriteReplacements(replacements, temp + position, from,
                                     length, escape_byte, &freqs);
  std::copy(temp, temp + result_length, from);
  delete [] temp;
  return result_length;
}

} // namespace bwtc
