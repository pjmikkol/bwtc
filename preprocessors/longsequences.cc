#include <cassert>
#include <cstdlib>

#include <algorithm>
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

  void Periods(std::vector<long_seq> *periods) {
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
    prev_hash_ = (((prev_hash_ - c_*data_[pos - 1]) << 8) +
                  data_[pos + block_length_ - 1]) % size_;
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
      if(StrEq(*prev, *curr)) {
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
      std::cout << "Frequency: " << it->count << " Length: "
                << it->length << " Position: " << it->position << "\n";
      for(uint64 i = it->position; i < it->position + it->length; ++i)
        std::cout << data_[i];
      std::cout << "\n-----------------------\n";
      ++s;
    }
    std::cout << "Total size of table: " << s;
    std::cout << "\n#########################\n";
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

  bool StrEq(const long_seq& s1, const long_seq& s2) {
    if (s1.length != s2.length) return false;
    for(int i = 0; i < s1.length; ++i) {
      if(data_[s1.position + i] != data_[s2.position + i])
        return false;
    }
    return true;
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
            /*(k == 0 && from[k] != from[start_pos]) ||
              (k != 0 && start_pos != k + block_length))*/
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
    std::cout << "Total matches in detecting long sequences: "
              << matches << "\n";
  seq_table.Periods(periods);
}

/* Replacement struct is used in lists which are inspected when writing *
 * replcements of sequences. */
struct replacement {
  explicit replacement(int64 loc, unsigned len, byte s) :
      location(loc), length(len), rank(s) {}

  int64 location;
  unsigned length;
  byte rank;
};

void DecreaseFreqs(FreqTable *freqs, unsigned *vals, unsigned times) {
  for(unsigned i = 0; i < 256; ++i) {
    if(vals[i]) {
      vals[i] *= times;
      freqs->Decrease(i, vals[i]);
    }
  }
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
    CalculateFrequencies(seq_freqs, data + s.position,
                         data + s.position + s.length);
    DecreaseFreqs(freqs, seq_freqs, s.count);
    repl[(data[s.position] << 8) + data[s.position + 1]]
        .push_back(replacement(s.position, s.length, j++));
  }
  if(j == 255) return 0xF000;
  std::partial_sort(periods->begin(), periods->begin() + (255-j), periods->end());
  unsigned i = 0;
  while(j < 255) {
    long_seq curr = (*periods)[i];
    if((*freqs)[j] >= curr.count*(curr.length - 1)) {
      break;
    } else {
      repl[(data[curr.position] << 8) + data[curr.position + 1]].
          push_back(replacement(curr.position, curr.length, j++));
      ++i;
    }
  }

  std::cout << "j=" << j << " i=" << i << "\n";

  //TODO: how to handle if only single replacement!!
  if(j-- <= 1) return 0xF000;
  if((*freqs)[j] > 0) {
    byte escape_index = j;
    while(j >= 0 && (*freqs)[j] > 0) {
      unsigned key = freqs->Key(j) << 8;
      for(unsigned k = 0; k < 256; ++k) {
        repl[key | k].push_back(replacement(-1,1,escape_index));
      }
      --j;
    }
    return (*freqs)[j];
  }
  return 0xF000;
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
  byte escape = DecideReplacements(&freqs, &periods, &long_seqs,
                                   replacements, from);
  bool escaping = escape < 0xff;

  return length;
}

} // namespace bwtc
