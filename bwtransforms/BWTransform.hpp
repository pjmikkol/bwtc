/**
 * @file BWTransform.hpp
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
 * Base classes (interfaces) for Burrows-Wheeler transform and its reversal.
 */

#ifndef BWTC_BWTRANSFORM_HPP_
#define BWTC_BWTRANSFORM_HPP_

#include <cassert>

#include <algorithm> // for reverse
#include <vector>

#include "../BWTBlock.hpp"
#include "../globaldefs.hpp"

namespace bwtc {

// TODO: delete old stuff
// If there is need for optimized memory management, then transformer
// needs to be connected to some manager-object  

/**
 * For implementing new algorithm for Burrows-Wheeler Transform one needs to
 * inherit BWTransform. After that BWTManager has to be modified.
 */
class BWTransform {
 public:
  BWTransform() {}
  virtual ~BWTransform() {}
  
  virtual
  void doTransform(byte *begin, uint32 length, std::vector<uint32>& LF) const = 0;

  virtual
  void doTransform(byte *begin, uint32 length, std::vector<uint32>& LF,
                   uint32 freqs[256]) const = 0;

  void doTransform(BWTBlock& block);
  void doTransform(BWTBlock& block, uint32 freqs[256]);

  virtual uint64 maxSizeInBytes(uint64 block_size) const = 0;
  virtual uint64 maxBlockSize(uint64 memory_budget) const = 0;
  virtual uint64 suggestedBlockSize(uint64 memory_budget) const = 0;

 private:
  BWTransform(const BWTransform&);
  const BWTransform& operator=(const BWTransform& );
};

// Block size and memory budget would probably be suitable parameters...
BWTransform* giveTransformer(char transform);

/**Computes the BWT from the suffix array and writes it to output.
 * 
 * Caller has to be sure that output has at least length of block_size + 1
 *
 * @param block Original string.
 * @param block_size Length of the original string.
 * @param suffix_array Siffx array of the block.
 * @param output BWT will be written to here.
 *
 * @return Position of the end-of-block character in output.
 */
template <typename Integer>
uint64 BWTFromSuffixArray(const byte* block, int64 block_size,
                          const Integer* suffix_array, byte* output)
{
  byte ch = '\0';
  Integer eob_position = 0;
  for (Integer rank = 0; rank < block_size + 1; ++rank) {
    Integer suffix = suffix_array[rank];
    assert(suffix <= block_size);
    if (suffix != 0) {
      ch = block[suffix - 1];
    } else {
      /* This is the end-of-block (eob) character, which cannot be
       * directly written since it has no code.
       * Instead, write a copy of the previous character here and
       * store this position in order to indicate its location. */
      eob_position = rank;
    }
    *output++ = ch;
  }
  return eob_position;
}

} //namespace bwtc

#endif
