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
 * Header for SequenceDetector-class.
 */

#ifndef BWTC_SEQUENCE_DETECTOR_H_
#define BWTC_SEQUENCE_DETECTOR_H_

#include <vector>
#include <list>
#include "../globaldefs.h"


namespace bwtc {
namespace long_sequences {


/**
 * Represents one interval of the input. Chunks are stored in list where
 * they are in the same order as they are in input. Besides that the
 * chunks representing the same string of bytes are linked together.
 * So we have multiple linked list in linked list.
 */
struct chunk {
  uint32 parent; /**<Link to link_table_ or long_table_ depending how long
                    the chunk is.*/
  uint32 position;
  std::list<chunk>::iterator prev;
  std::list<chunk>::iterator next;
};

/**
 * Represents single entry in link_table_
 */
struct link_struct {
  /**<Iterator to chunk_list_. Length of this entry
     is the size of window at most. */
  std::list<chunk>::iterator short_chunk; 
  
  //TODO: is there need to save memory in this. length is always small and
  //      needs at most few bytes. 
  /**<Length of the chunk. Always lesser or equal to window
     size. Smaller when we are dealing with periods. */
  uint32 length; 

  /**<Next entry in link_table_ which has the same hash value. */
  int next;

  /**<Count of this chunk idetified so far. */  
  uint32 count;
  
  /**<If some entry of this chunk is later found to
     be a part of some long sequence then this integer
     is the index of that long sequence in long_table_. */
  int related_long_seq; 
  
  /**<Offset position for this chunk in the
     related_long_seq.
  */
  int place_in_long; 
};

/**
 * Represents single entry in SequenceDetector<Hasher>#long_table_.
 */
struct long_seq {
  long_seq(std::list<chunk>::iterator long_chunk, uint32 l, uint32 c) :
      head(long_chunk), length(l), count(c) {}
  
  bool operator<(const long_seq& ll) const {
    return count*length < ll.count*ll.length;
  }
  
  std::list<chunk>::iterator head;
  uint32 length;
  uint32 count;
};


/**
 * SequenceDetector finds repeating sequences from the given input.
 *
 * Roots of the algorithm lies on the Karp-Rabin-string search algorithm and
 * Bentley-McIlroy compression algorithm. 
 */
template <typename Hasher>
class SequenceDetector {
  static const int kMinPeriod = 16;
  static const int kWindowSize = 64;



  
 public:
  // note that only 32-bit lengths are supported
  // period_threshold (kMinPeriod) and window_size are compile time constants
  SequenceDetector(byte *from, uint64 length, byte *freqs) {
    
  }
  
 private:
  Hasher h_; /**<Hasher-object which computes the rolling-hash function.
                @see hash_functions */
  uint32 *h_table_; /**<Table indexed by hash-values. Contains an index of
                       link-table of the chain of this hash-value.*/
  uint32 h_table_size_; //TARVITAANKO?
  std::vector<link_struct> link_table; /**<Link-table which holds the chains
                                          of different hash-values. Each member
                                          points to the head of the
                                          sequence-list. */
  std::vector<long_seq> long_table_; /**<Analogy of link-table for long
                                        sequences. */
  std::list<chunk> chunk_list_; /**<Members of this list represent one chunk
                                   of input data. Each chunk has correspondent
                                   entry in either link_table or long_table_.*/


};

} //namespace long_sequences
} //namespace bwtc

#endif
