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

#ifndef BWTC_BLOCK_H_
#define BWTC_BLOCK_H_

#include <cassert>

#include <vector>

#include "globaldefs.h"

namespace bwtc {

/**
 * MainBlock models one block of a data which is en-/decoded as an unit.
 *
 * On a compressor pipeline the preprocessor produces MainBlocks which are
 * then transformed into a sequence of context blocks by BWTransform. 
 * MainBlock has two  arrays: {@link #block_} and {@link #stats_}. Both are
 * allocated and deallocated by the client. BlockManager-class exists for this
 * task.
 *
 * @see BlockManager
 * @see BWTransform
 * @see PreProcessor
 */
class MainBlock {
 public:
  /**
   * Creates MainBlock-object from given allocated resources.
   *
   * @param block Pointer to vector which will be used as an storage
   *              for the data the MainBlock-object holds ({@link #block_}).
   * @param stats Pointer to vector which will be used to hold frequencies
   *              of the characters ({@link #stats_}).
   * @param filled Range [0,filled) of the block contains meaningful data .
   */
  MainBlock(std::vector<byte>* block, std::vector<uint64>* stats,
            uint64 filled);
  /**
   * Destructor which doesn't free any resources the object has used!
   */
  ~MainBlock();

  /**
   * Size of actual data in object.
   *
   * @return number of bytes of actual data in object
   */
  inline uint64 Size() { return filled_; }

  /**
   * Start of the data of the object.
   *
   * @return memory address for the first byte of the data the object holds
   */
  inline byte* begin()  { return &(*block_)[0]; }

  /**
   * End of the data of the object.
   *
   * @return memory address one past the last byte of the data the object holds
   */
  inline byte* end() { return &(*block_)[filled_]; }

  /**
   * Appends single byte to the end of filled area.
   *
   * @param toAppend byte which will be appended
   */
  inline void Append(byte toAppend) {
    assert(filled_ < block_->size());
    (*block_)[filled_++] = toAppend;
  }
  
  std::vector<byte> *block_; /**<Contains the actual data of MainBlock-object*/
  std::vector<uint64> *stats_;/**<Contains frequencies for each byte in block_.
                                 After Burrows-Wheeler transform this holds the
                                 sizes of context blocks of a single byte in
                                 the same order as they appear in transformed
                                 data.*/
  uint64 filled_;/**<Range [0, filled_) in block_ will hold the relevant data.*/


 private:
  MainBlock& operator=(const MainBlock& b);
  MainBlock(const MainBlock&);

};

} // namespace bwtc

#endif
