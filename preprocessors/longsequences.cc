/**
 * @file longsequences.cc
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
 * Implementation for finding repetitive sequences and replacing them.
 */
#include <cassert>
#include <cstdlib>

#include <algorithm>
#include <list>
#include <utility>
#include <vector>

#include "longsequences.h"
#include "preprocessor.h"
//#include "sequence_detector.h"
#include "../coders.h"
#include "../globaldefs.h"
#include "../utils.h"

#include <iostream>

namespace bwtc {

namespace long_sequences {


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
 * replacements of sequences. */
struct replacement {
  explicit replacement(int64 loc, unsigned len, byte s) :
      location(loc), length(len), rank(s) {}

  int64 location;
  unsigned length;
  byte rank; /* ????? */
};

/* We do the comparison of the srings backwards, because it seems to be *
 * faster the forwards. */
bool SeqEq(const long_seq& s1, const long_seq& s2, byte *data) {
  if (s1.length != s2.length) return false;
  byte *f = data + s1.position + s1.length - 1;
  byte *g = data + s2.position + s1.length - 1;
  const byte *start = data + s1.position - 1;
  while(f != start) if(*f-- != *g--) return false;
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

void DetectSequences(byte *from, uint64 length, int memory_constraint,
                     unsigned block_length, int threshold, uint64 *freqs,
                     LongSequences *long_seqs, std::vector<long_seq> *periods)
{

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

}
/**
 * Replaces sequences in given byte-array,
 *
 * @param rpls Array of size 65536 which holds the sequences to be replaced.
 *             Each sequence is stored in list which the table holds. Lists
 *             are indexed by the first two bytes of the replaceable sequence.
 * @param to Byte-array where the result is written.
 * @param from Byte-array of the unmodified input string.
 * @param length Length of the input.
 * @param escape_byte Escape byte used in writing the replacements.
 * @param freqs Pointer to frequency table, which holds the frequencies of the
 *              characters in input.
 * @return Bytes used for writing the encoded message.
 */
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
                         unsigned window_size, int threshold)
{
  /*  using namespace long_sequences;
  assert(length > 0);
  assert(from);
  assert(length > window_size);
  assert(memory_constraint > 1);

  uint64 frequencies[256] = {0};
  
  DetectSequences(from, length, memory_constraint, window_size, threshold,
                  frequencies, &long_seqs, &periods);
  FreqTable freqs(frequencies);
  std::list<replacement> replacements[65536];
  unsigned escape_index = DecideReplacements(&freqs, &periods, &long_seqs,
                                             replacements, from);
  bool escaping = escape_index < 255;
  byte escape_byte = escaping ? freqs.Key(escape_index) : 0;
  byte *temp = new byte[2*length];
  unsigned position = 0;*/
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
  /*  unsigned prev_s = 0xF000;
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
  return result_length;*/
}

} // namespace bwtc
