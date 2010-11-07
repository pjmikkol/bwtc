/**
 * @file sequence_detector-inl.h
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
 * Template based implementation of the SequenceDetector-class.
 */

#ifndef BWTC_SEQUENCE_DETECTOR_INL_H_
#define BWTC_SEQUENCE_DETECTOR_INL_H_

#include <algorithm>
#include <cassert>
#include <vector>
#include "../globaldefs.h"
#include "longsequences.h"
#include "sequence_detector.h"
#include <utility>

namespace bwtc {
namespace long_sequences {

/**
 * SequenceDetector finds repeating sequences from the given input.
 *
 * Roots of the algorithm lies on the Karp-Rabin-string search algorithm and
 * Bentley-McIlroy compression algorithm. 
 */
template <typename Hasher>
class SequenceDetector {

 public:
  // note that only 32-bit lengths are supported
  // period_threshold (kMinPeriod) and window_size are compile time constants
  SequenceDetector(byte *from, uint32 h_table_size, uint32 *freqs,
                   std::vector<chunk> *chunks, std::vector<uint32> *h_table,
                   uint32 win_size) :
      chunks_(*chunks), h_table_(*h_table), window_size(win_size), source_(from)
  {
    h_.Initialize(h_table_size, window_size);
    h_table_.resize(h_.Size());
    std::fill(h_table_.begin(), h_table_.end(), 0);
    chunks_.reserve(h_table_size);
    position_ = 0;
    freqs_ = freqs;
  }

  /**
   *
   */
  uint32 Count(const byte *string, uint32 len) {
    return h_table_[h_.InitValue(string, len)];
  }

  /**
   * Scans the next length bytes of the given source.
   *
   * @param length how many bytes we proceed
   * @return how many bytes we proceeded.
   */
  uint32 ScanAndStore(uint32 length) {
    assert(length >= window_size);
    length += position_;
    int64 h_val = h_.InitValue(source_ + position_, window_size);
    InsertToHashTable(h_val, position_);
    CalculateFrequencies(freqs_, source_ + position_,
                         source_ + position_ + window_size);
    position_ += window_size;
    while(position_ <= length - 2*window_size - 1) {
      uint32 old_pos = position_;
      ScanAndCompare(window_size);
      CalculateFrequencies(freqs_, source_ + old_pos, source_ + position_);
    }
    CalculateFrequencies(freqs_, source_ + position_, source_ + length);
    if(position_ <= length - window_size - 1) {
      ScanAndCompare((length - position_) - window_size + 1);
    }
    return position_; //TODO: fix this
  }

  size_t ChunksCount() {
    return chunks_.size();
  }

 private:
  /**
   *
   */
  void InsertToHashTable(uint32 h_val, uint32 pos) {
    ++h_table_[h_val];
    chunks_.push_back(chunk(pos, h_val));
  }
  void InsertToHashTable(chunk c) {
    ++h_table_[c.hash_value];
    chunks_.push_back(c);
  }

  /**
   *
   */
  void ScanAndCompare(uint32 b) {
    // TODO: memory allocation optimization if needed
    std::vector<chunk> temp_chunks;
    temp_chunks.reserve(b);
    int64 h_val = h_.InitValue(source_ + position_, window_size);
    temp_chunks.push_back(chunk(position_, h_val));
    for(uint32 i = 1; i < b; ++i) {
      h_val = h_.Update(source_[position_], source_[position_ + window_size]);
      temp_chunks.push_back(chunk(++position_, h_val));
    }
    std::sort(temp_chunks.begin(), temp_chunks.end(), cmp_chunk(&h_table_[0]));
    std::pair<uint32, uint32> m = MaxDuplicates(temp_chunks);
    for(size_t i = m.first; i < m.first + m.second; ++i) {
      InsertToHashTable(temp_chunks[i]);
    }
    position_ = temp_chunks[m.first + m.second-1].position + window_size;
  }

  /**
   *
   *
   * @return Pair where second element is the count of duplicates and
   *         the first element tells the first index of the duplicate values.
   */
  std::pair<uint32, uint32> MaxDuplicates(const std::vector<chunk>& chunks) {
    uint32 ind = chunks.size() - 1, count = 1, max_pos = ind;
    for(int i =  ind; i > 0; --i) {
      uint32 d_count = 1;
      while(i > 0 && chunks[i].hash_value == chunks[i-1].hash_value) {
        ++d_count; --i;
      }
      if(d_count > count) { count = d_count; ind = i; }
      else if(h_table_[chunks[max_pos].hash_value] == 0 &&
              chunks[max_pos].position < chunks[i].position) max_pos = i;
    }
    if(count == 1) ind = max_pos;
    return std::make_pair(ind, count);
  }


  /**<Hasher-object which computes the rolling-hash function.
     @see hash_functions */
  Hasher h_; 
  //TODO: Used at least in first phase ..
  std::vector<chunk>& chunks_;
  /**<Table indexed by hash-values. Contains an index of
     link-table of the chain of this hash-value.*/
  std::vector<uint32>& h_table_;
  /**<Pointer to the current spot of the sequence.*/
  byte *source_;
  /**< */
  uint32 position_;
  uint32 *freqs_;
  uint32 window_size;
};

} //namespace long_sequences
} //namespace bwtc

#endif
