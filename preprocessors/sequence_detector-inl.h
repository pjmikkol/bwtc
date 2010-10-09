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
#include <vector>
#include "../globaldefs.h"
#include "longsequences.h"
#include "sequence_detector.h"
#include <utility>

namespace bwtc {
namespace long_sequences {

/**
 * Namespace for struct-definitions and other helpers used in SequenceDetector.
 */
namespace {

//TODO: name...
struct chunk {
  chunk(uint32 p, uint32 h) : position(p), hash_value(h) {}
  uint32 position;
  uint32 hash_value;
};

struct cmp_chunk {
  cmp_chunk(uint32 *c) : counts(c) {}
  bool operator()(const chunk& a, const chunk& b) {
    if (counts[a.hash_value] < counts[b.hash_value])
      return true;
    else if (counts[a.hash_value] > counts[b.hash_value])
      return false;
    else if (a.hash_value < b.hash_value)
      return true;
    else if (a.hash_value > b.hash_value)
      return false;
    else return a.position < b.position;
  }
  uint32 *counts;
};

} //empty namespace

#ifndef BWTC_SEQ_DET_SCAN_CONSTS_
const uint32 kMinPeriod = 16;
const uint32 kWindowSize = 32;
#endif

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
  SequenceDetector(byte *from, uint32 h_table_size, uint32 *freqs) {
    source_ = from;
    h_.Initialize(h_table_size, kWindowSize);
    h_table_ = new uint32[h_.Size()];
    std::fill(h_table_, h_table_ + h_.Size(), 0);
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
    assert(length >= kWindowSize);
    length += position_;
    int64 h_val = h_.InitValue(source_ + position_, kWindowSize);
    InsertToHashTable(h_val, position_);
    CalculateFrequencies(freqs_, source_ + position_,
                         source_ + position_ + kWindowSize);
    position_ += kWindowSize;
    while(position_ <= length - 2*kWindowSize - 1) {
      uint32 old_pos = position_;
      ScanAndCompare(kWindowSize);
      CalculateFrequencies(freqs_, source_ + old_pos, source_ + position_);
    }
    CalculateFrequencies(freqs_, source_ + position_, source_ + length);
    if(position_ <= length - kWindowSize - 1) {
      ScanAndCompare((length - position_) - kWindowSize + 1);
    }
    return position_; //TODO: fix this
  }

  ~SequenceDetector() {
    delete [] h_table_;
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
    int64 h_val = h_.InitValue(source_ + position_, kWindowSize);
    temp_chunks.push_back(chunk(position_, h_val));
    for(uint32 i = 1; i < b; ++i) {
      h_val = h_.Update(source_[position_], source_[position_ + kWindowSize]);
      temp_chunks.push_back(chunk(++position_, h_val));
    }
    std::sort(temp_chunks.begin(), temp_chunks.end(), cmp_chunk(h_table_));
    std::pair<uint32, uint32> m = MaxDuplicates(temp_chunks);
    for(size_t i = m.first; i < m.first + m.second; ++i) {
      InsertToHashTable(temp_chunks[i]);
    }
    position_ = temp_chunks[m.first + m.second-1].position + kWindowSize;
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
  /**<Table indexed by hash-values. Contains an index of
     link-table of the chain of this hash-value.*/
  uint32 *h_table_;
  /**<Pointer to the current spot of the sequence.*/
  byte *source_;
  //TODO: Used at least in first phase ..
  std::vector<chunk> chunks_;
  /**< */
  uint32 position_;
  uint32 *freqs_;
  

  
  uint32 h_table_size_; //TARVITAANKO?
};

} //namespace long_sequences
} //namespace bwtc

#endif
