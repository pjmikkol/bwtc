/**
 * @file BWTransform.cpp
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
 * Implementation for the base class of Burrows-Wheeler transformers.
 */

#include <cassert>
#include <algorithm>
#include <vector>

#include "../globaldefs.hpp"
#include "../BWTBlock.hpp"
#include "BWTransform.hpp"
#include "SA-IS-bwt.hpp"
#include "Divsufsorter.hpp"

namespace bwtc {

void BWTransform::doTransform(BWTBlock& block) {
  std::reverse(block.begin(), block.end());
  byte next = *block.end();
  *block.end() = 0;

  doTransform(block.begin(), block.size() + 1, block.LFpowers());
  block.setTransformed(true);

  *(block.begin() + block.LFpowers()[0]) = *block.end();
  
  *block.end() = next;
}

void BWTransform::doTransform(BWTBlock& block, uint32 freqs[256]) {
  std::reverse(block.begin(), block.end());
  byte next = *block.end();
  *block.end() = 0;

  doTransform(block.begin(), block.size() + 1, block.LFpowers(), freqs);
  block.setTransformed(true);

  *(block.begin() + block.LFpowers()[0]) = *block.end();

  
  *block.end() = next;
}

BWTransform* giveTransformer(char transform) {
  (void) transform;
  if(transform != 's') {
    if(verbosity > 1) std::clog << "Using divsufsort for calculating BWT.\n";
    return new Divsufsorter();
  }
  else {
    if(verbosity > 1) std::clog << "Using sais for calculating BWT.\n";
    return new SAISBWTransform();
  }
}

}
