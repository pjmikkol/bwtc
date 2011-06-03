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

#include <cassert>

#include <vector>

#include "bw_transform.h"
#include "sa-is-bwt.h"
#include "sais.hxx"
#include "../MainBlock.hpp"
#include "../globaldefs.hpp"

namespace bwtc {

SAISBWTransform::SAISBWTransform() {}

std::vector<byte>* SAISBWTransform::DoTransform(uint64 *eob_byte) {
  if(!current_block_) return NULL;

  int block_size = current_block_->size();
  current_block_->append(0);
  byte *block = current_block_->begin();
  /* Whole transformation is done in single pass. */
  current_block_ = NULL;
  std::vector<byte> *result = AllocateMemory(block_size + 1);
  std::vector<int> suffix_array(block_size + 1);
  *eob_byte = saisxx_bwt(block, &(*result)[0], &suffix_array[0], block_size + 1);
  /* Slower implementation of the same algorithm (unsigned are supported)
  sa_is::SA_IS_ZeroInclude(block, &suffix_array[0], block_size + 1, 256);
    *eob_byte = BwtFromSuffixArray(block, block_size, &suffix_array[0],
                               &(*result)[0]); */
  return result;
}

} //namespace bwtc
