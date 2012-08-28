/**
 * @file sa-is-bwt.hpp
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
 * Header for SA-IS BWT-algorithm.
 */

#ifndef BWTC_SAIS_BWT_HPP_
#define BWTC_SAIS_BWT_HPP_

#include <cassert>

#include <algorithm>
#include <iostream>
#include <vector>

#include "BWTransform.hpp"
#include "../globaldefs.hpp"

using bwtc::uint64;
using bwtc::int64;
using bwtc::uint32;
using bwtc::byte;
using bwtc::Max;

namespace bwtc {

/**
 * Burrows-Wheeler transform using the SA-IS algorithm desribed in
 * "Two Efficient Algorithms for Linear Suffix Array Construction"
 * by Nong, Zhang & Chan
 */
class SAISBWTransform : public BWTransform {
 public:
  SAISBWTransform();
  virtual ~SAISBWTransform() {}
  void
  doTransform(byte *begin, uint32 length, std::vector<uint32>& LFpowers) const;

  void
  doTransform(byte *begin, uint32 length, std::vector<uint32>& LFpowers,
              uint32* freqs) const;

  /* The following values aren't correct */
  virtual uint64 maxSizeInBytes(uint64) const { return 0; }
  virtual uint64 maxBlockSize(uint64) const { return 0; }
  virtual uint64 suggestedBlockSize(uint64) const { return 0; }
  
};
} // namespace bwtc


#endif
