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
#include "../MainBlock.hpp"
#include "../globaldefs.hpp"
#include "../Profiling.hpp"

namespace bwtc {

SAISBWTransform::SAISBWTransform() {}

std::vector<byte>* SAISBWTransform::doTransform(std::vector<uint32>& LFpowers) {
  PROFILE("SAISBWTransform::doTransform");
  if(!m_currentBlock) return 0;

  int block_size = m_currentBlock->size();
  m_currentBlock->append(0);
  byte *block = m_currentBlock->begin();
  /* Whole transformation is done in single pass. */
  m_currentBlock = 0;
  std::vector<byte> *result = allocateMemory(block_size + 1);
  std::vector<int> suffix_array(block_size + 1);
  saisxx_bwt(block, &(*result)[0], &suffix_array[0], block_size + 1, LFpowers);
  /* Slower implementation of the same algorithm (unsigned are supported)
  sa_is::SA_IS_ZeroInclude(block, &suffix_array[0], block_size + 1, 256);
    *eob_byte = BwtFromSuffixArray(block, block_size, &suffix_array[0],
                               &(*result)[0]); */
  return result;
}

} //namespace bwtc
