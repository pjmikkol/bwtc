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

#include "../globaldefs.h"
#include "bw_transform.h"
#include "dcbwt.h"
#include "sa-is-bwt.h"

namespace bwtc {

/* We assume that the context-block of sentinel char is at the front of *
 * transform.                                                           */
//TODO: Allow to take array from preprocessor to parameter since they have
//      already built stats
void BWTransform::BuildStats() {
  std::fill(current_block_->m_stats->begin(), current_block_->m_stats->end(), 0);
  (*current_block_->m_stats)[0] = 1;
  for(uint64 i = 0; i < current_block_->size(); ++i)
    (*current_block_->m_stats)[(*current_block_->m_block)[i] + 1]++;
}

std::vector<byte>* BWTransform::AllocateMemory(uint64 size) {
  return new std::vector<byte>(size);
}

BWTransform* GiveTransformer(char transform) {
  /*When there are multiple ways to do transform this is to place to add them*/
  if(transform == 'd')
    return new DCBWTransform(8);
  else
    return new SAISBWTransform();
}

}
