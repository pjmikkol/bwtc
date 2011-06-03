/**
 * @file MainBlock.cpp
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
 * Implementation for the MainBlock-class.
 */

#include <vector>
#include <boost/cstdint.hpp>

#include "MainBlock.hpp"
#include "globaldefs.hpp"

namespace bwtc {

MainBlock::MainBlock(std::vector<byte>* block, std::vector<uint64>* stats,
                     uint64 filled) : 
    m_block(block), m_stats(stats), m_filled(filled) {}

MainBlock::~MainBlock() {}

} //namespace bwtc
