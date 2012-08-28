/**
 * @file Divsufsorter.hpp
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
 * Header for divsufsort-algorithm (implemented by Yuta Mori). 
 */

#ifndef BWTC_DIVSUFSORT_BWT_HPP_
#define BWTC_DIVSUFSORT_BWT_HPP_

#include <cassert>

#include <algorithm>
#include <vector>

#include "BWTransform.hpp"
#include "../globaldefs.hpp"
#include "../Profiling.hpp"

#include "divsufsort.h"

using bwtc::uint64;
using bwtc::int64;
using bwtc::uint32;
using bwtc::byte;
using bwtc::Max;

namespace bwtc {

class Divsufsorter : public BWTransform {
 public:
  Divsufsorter() {}
  virtual ~Divsufsorter() {}

  void
  doTransform(byte *begin, uint32 length, std::vector<uint32>& LFpowers) const {
    PROFILE("Divsufsorter::doTransform");
    divbwt(begin, begin, 0, length, &LFpowers[0], LFpowers.size());
  }

  void
  doTransform(byte *begin, uint32 length, std::vector<uint32>& LFpowers,
              uint32 *freqs) const {
    PROFILE("Divsufsorter::doTransform");
    divbwtf(begin, begin, 0, length, &LFpowers[0], LFpowers.size(), freqs);
  }

  /* The following values aren't correct */
  virtual uint64 maxSizeInBytes(uint64) const { return 0; }
  virtual uint64 maxBlockSize(uint64) const { return 0; }
  virtual uint64 suggestedBlockSize(uint64) const { return 0; }
  
};
} // namespace bwtc


#endif
