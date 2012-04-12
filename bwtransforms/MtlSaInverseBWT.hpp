/**
 * @file MtlSaInverseBWT.hpp
 * @author Dominik Kempa <dominik.kempa@cs.helsinki.fi>
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
 * Header for MTL-SA algorithm for inverting BWT.
 */

#ifndef MTL_SA_INVERSE_BWT_HPP_
#define MTL_SA_INVERSE_BWT_HPP_

#include <vector>

#include "../globaldefs.hpp"
#include "InverseBWT.hpp"

namespace bwtc {

/**
 * Inverse Burrows-Wheeler transform using the MTL-SA algorithm described in
 * "Slashing the Time for BWT Inversion" by Karkkainen, Kempa and Puglisi.
 */
class MtlSaInverseBWTransform : public InverseBWTransform {
 public:
  MtlSaInverseBWTransform() {}
  virtual ~MtlSaInverseBWTransform() {}
  virtual uint64 maxBlockSize(uint64 memory_budget) const;
  virtual std::vector<byte>* doTransform(const byte* source_bwt,
                                         uint64 bwt_size,
                                         const std::vector<uint32> &LFpowers);
};

} //namespace bwtc

#endif
