/**
 * @file BWTransform.cpp
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
 * Implementation for the base class of Burrows-Wheeler transformers.
 */

#include <cassert>

#include <vector>

#include "../globaldefs.hpp"
#include "BWTransform.hpp"
#include "SA-IS-bwt.hpp"

namespace bwtc {

/* We assume that the context-block of sentinel char is at the front of
 * transform i.e. sentinel char is the smallest character in alphabet.
 */
//TODO: Allow to take array from preprocessor to parameter since they have
//      already built stats
void BWTransform::buildStats() {
  std::fill(m_currentBlock->m_stats->begin(), m_currentBlock->m_stats->end(), 0);
  (*m_currentBlock->m_stats)[0] = 1;
  for(uint64 i = 0; i < m_currentBlock->size(); ++i)
    (*m_currentBlock->m_stats)[(*m_currentBlock->m_block)[i] + 1]++;
}

std::vector<byte>* BWTransform::allocateMemory(uint64 size) {
  return new std::vector<byte>(size);
}

BWTransform* giveTransformer(char transform) {
    return new SAISBWTransform();
}

}
