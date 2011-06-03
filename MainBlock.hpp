/**
 * @file MainBlock.hpp
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
 * Header for MainBlock-class.
 */

#ifndef BWTC_MAINBLOCK_HPP_
#define BWTC_MAINBLOCK_HPP_

#include <cassert>
#include <vector>

#include "globaldefs.hpp"

namespace bwtc {

/**
 * MainBlock models one block of a data which is en-/decoded as an unit.
 *
 * On a compressor pipeline the preprocessor produces MainBlocks which are
 * then transformed into a sequence of context blocks by BWTransform. 
 * MainBlock has two  arrays: {@link #m_block} and {@link #m_stats}. Both are
 * allocated and deallocated by the user. BlockManager-class exists for this
 * task. It is possible that the MainBlock is greater than data it holds.
 *
 * @see BlockManager
 * @see BWTransform
 * @see PreProcessor
 */
class MainBlock {
 public:
  /**
   * Creates MainBlock-object from given already allocated resources.
   *
   * @param block Pointer to vector which will be used as an storage
   *              for the data the MainBlock-object holds ({@link #m_block}).
   * @param stats Pointer to vector which will be used to hold frequencies
   *              of the characters ({@link #m_stats}).
   * @param filled Range [0,filled) of the block contains meaningful data.
   */
  MainBlock(std::vector<byte>* block, std::vector<uint64>* stats,
            uint64 filled);
  /**
   * Destructor. Doesn't free any resources the object has used!
   */
  ~MainBlock();

  /**
   * Size of actual data in object.
   *
   * @return number of bytes of actual data in object
   */
  inline uint64 size() { return m_filled; }

  /**
   * Start of the data of the object.
   *
   * @return memory address for the first byte of the data the object holds
   */
  inline byte* begin()  { return &(*m_block)[0]; }

  /**
   * End of the data of the object.
   *
   * @return memory address one past the last byte of the data the object holds
   */
  inline byte* end() { return &(*m_block)[m_filled]; }

  /**
   * Appends single byte to the end of filled area.
   *
   * @param toAppend byte which will be appended
   */
  inline void append(byte toAppend) {
    assert(m_filled < m_block->size());
    (*m_block)[m_filled++] = toAppend;
  }
  
  /**<Contains the actual data of MainBlock-object*/
  std::vector<byte> *m_block; 
  /**<Contains frequencies for each byte in m_block. After Burrows-Wheeler
      transform this holds the sizes of context blocks of a single byte in
      the same order as they appear in transformed data.*/
  std::vector<uint64> *m_stats;
  /**<Range [0, m_filled) in m_block will hold the relevant data.*/
  uint64 m_filled;


 private:
  MainBlock& operator=(const MainBlock& b);
  MainBlock(const MainBlock&);

};

} // namespace bwtc

#endif
