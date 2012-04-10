/**
 * @file MTL-SA-ibwt.hpp
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

#ifndef BWTC_MTL_SA_IBWT_HPP_
#define BWTC_MTL_SA_IBWT_HPP_

#include <vector>

#include "../globaldefs.hpp"
#include "InverseBWT.hpp"

namespace bwtc {

/**
 * Inverse Burrows-Wheeler transform using the MTL-SA algorithm described in
 * "Slashing the Time for BWT Inversion" by Karkkainen, Kempa and Puglisi.
 */
class MTL_SA_InverseBWTransform : public InverseBWTransform {
 public:
  MTL_SA_InverseBWTransform() {}
  virtual ~MTL_SA_InverseBWTransform() {}
  virtual uint64 maxBlockSize(uint64 memory_budget) const;
  virtual std::vector<byte>* doTransform(const byte* source_bwt,
                                         uint64 bwt_size,
                                         uint64 eob_position);
 private:
 /** Computes the number of occurrences of each possible pair of bytes in the 
  *  original string. Additionally, saves these pairs (in the BWT order) into
  *  'data' array (according to predetermined layout, see doTransform for
  *  details).
  *
  *  @param[in] bwt           Burrows-Wheeler transformed string.
  *  @param[in] bwt_size      Length of bwt.
  *  @param[in] eob_position  Position of EOB character in bwt.
  *  @param[out] data         An array holding the character pairs.
  *  @param[out] pairs_count  Character pairs counts.
  *  @return                  Indices of rows of the BWT matrix that end with
  *                           a pair containing the EOB symbol.
  */
  std::pair<uint32, uint32> computePairs(const byte* bwt, uint64 bwt_size,
                                        uint64 eob_position,
                                        std::vector<uint32> &data,
                                        std::vector<uint32> &pairs_count);

 /** Compute the LF^2[i] for each position i in the BWT and save to 'data'.
  *
  *  @param[in] bwt_size      Length of bwt.
  *  @param[in] EOB_pairs     Indices of rows of the BWT matrix that end with
  *                           a pair containing the EOB symbol.
  *  @param[in] pairs_count   Number of occurrences in the original string of
  *                           each pair.
  *  @param[in,out] data      An array with character pairs storing also the
  *                           output of this function.
  */
  void computeLF2(uint64 bwt_size, std::pair<uint32, uint32> &EOB_pairs,
                  std::vector<uint32> &pairs_count, std::vector<uint32> &data);
};

} //namespace bwtc

#endif
