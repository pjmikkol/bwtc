/**
 * @file sa-is-bwt.cpp
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
 * Wrapper for the SA-IS implementation.
 */

#include <cassert>

#include <vector>

#include "BWTransform.hpp"
#include "SA-IS-bwt.hpp"
#include "sais.hxx"
#include "../globaldefs.hpp"
#include "../Profiling.hpp"

namespace bwtc {

SAISBWTransform::SAISBWTransform() {}

void SAISBWTransform::
doTransform(byte *begin, uint32 length, std::vector<uint32>& LFpowers) const {
  PROFILE("SAISBWTransform::doTransform");
  std::vector<int> SA(length);
  saisxx_bwt(begin, begin, &SA[0], (int)length, LFpowers);
}

void SAISBWTransform::
doTransform(byte *begin, uint32 length, std::vector<uint32>& LFpowers,
            uint32 *freqs) const {
  PROFILE("SAISBWTransform::doTransform");
  std::vector<int> SA(length);
  saisxx_bwt(begin, begin, &SA[0], (int)length, LFpowers, 256, freqs);
}

} //namespace bwtc
