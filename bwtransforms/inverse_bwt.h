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

#ifndef BWTC_INVERSE_BWT_H_
#define BWTC_INVERSE_BWT_H_

#include <vector>

#include "../globaldefs.hpp"

namespace bwtc {
/**************************************************************************
 * Base class for Inverse Burrows-Wheeler transforms. This class has some *
 * common functionality with BWTransform at the moment (especially memory *
 * allocation mechanism) so merging these two together is an option.      *
 *                                                                        *
 * At the moment biggest differece is that in inverse transform we aren't *
 * prepared to make the transform in pieces, so Connecting to the block   *
 * before transform is not required as it is in forward-transform.        *
 **************************************************************************/

class InverseBWTransform {
 public:
  InverseBWTransform() {}
  virtual ~InverseBWTransform() {}
  virtual uint64 MaxBlockSize(uint64 memory_budget) const = 0;
  virtual std::vector<byte>* DoTransform(const byte* source_bwt,
                                         uint64 bwt_size,
                                         uint64 eob_position) = 0;

 protected:
  /* This is here for the sake of memory management */
  virtual std::vector<byte>* AllocateMemory(uint64 block_size);
};

/**************************************************************************
 * FastInverseBWTransformer: This is the algorithm mergeTL in Sewards     *
 * paper with some additional heuristics to be able to handle very large  *
 * blocks.                                                                *
 **************************************************************************/
class FastInverseBWTransform : public InverseBWTransform {
 public:
  FastInverseBWTransform() {}
  virtual ~FastInverseBWTransform() {}
  virtual uint64 MaxBlockSize(uint64 memory_budget) const;
  virtual std::vector<byte>* DoTransform(const byte* source_bwt,
                                         uint64 bwt_size, uint64 eob_position);
 private:
  static const int64 kMemoryOverhead = 1 << 20;
};

/* For example memory budget would be a good parameter.. */
InverseBWTransform* GiveInverseTransformer();

} //namespace bwtc
#endif
