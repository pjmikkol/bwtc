/**
 * @file InverseBWT.hpp
 * @author Pekka Mikkola <pmikkol@gmail.com>
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
 * Header for inverse of BWT.
 */

#ifndef BWTC_INVERSE_BWT_HPP_
#define BWTC_INVERSE_BWT_HPP_

#include <vector>

#include "../globaldefs.hpp"
#include "../BWTBlock.hpp"

namespace bwtc {
/**
 * Base class for Inverse Burrows-Wheeler transforms. This class has some
 * common functionality with BWTransform at the moment (especially memory
 * allocation mechanism) so merging these two together is an option.
 *
 * At the moment biggest differece is that in inverse transform we aren't
 * prepared to make the transform in pieces, so Connecting to the block
 * before transform is not required as it is in forward-transform.
 */
class InverseBWTransform {
 public:
  virtual ~InverseBWTransform() {}
  virtual uint64 maxBlockSize(uint64 memory_budget) const = 0;

  virtual void doTransform(byte *bwt, uint32 n,
                           const std::vector<uint32>& LFpow) = 0;

  void doTransform(BWTBlock& block);

};

/**
 * This is the algorithm mergeTL in Sewards' paper with some additional
 * heuristics to be able to handle very large blocks.
 *
 * For original implementation see 
 * @see http://code.google.com/p/dcs-bwt-compressor/ 
 */
class FastInverseBWTransform : public InverseBWTransform {
 public:
  FastInverseBWTransform() {}
  virtual ~FastInverseBWTransform() {}
  virtual uint64 maxBlockSize(uint64 memory_budget) const;
  virtual void doTransform(byte* source_bwt,
                           uint32 bwt_size,
                           const std::vector<uint32>& LFpowers);
 private:
  static const int64 kMemoryOverhead = 1 << 20;
};

/* For example memory budget would be a good parameter.. */
InverseBWTransform* giveInverseTransformer();

} //namespace bwtc
#endif
