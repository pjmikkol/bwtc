/**
 * @file BWTBlock.hpp
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
 * Header for BWT-block. BWT-block is part of preprocessed data, that will
 * be transformed in single block. Each BWT-block will be encoded independently.
 */

#ifndef BWTC_BWTBLOCK_HPP_
#define BWTC_BWTBLOCK_HPP_

#include "globaldefs.hpp"

#include <cassert>
#include <vector>

namespace bwtc {

class BWTBlock {
 public:
  BWTBlock(byte *data, uint32 length, bool isTransformed);
  BWTBlock(const BWTBlock& b);
  BWTBlock& operator=(const BWTBlock& b);

  /**This should only be used by BWTransform-/InverseBWT-classes. */
  void setTransformed(bool transformed) {
    assert(m_isTransformed != transformed);
    m_isTransformed = transformed;
  }

  bool isTransformed() const { return m_isTransformed; }
  size_t size() const { return m_length; }
  byte* begin() { return m_begin; }
  const byte* begin() const { return m_begin; }
  byte* end() { return m_begin + m_length; }
  const byte* end() const { return m_begin + m_length; }
  std::vector<uint32>& LFpowers() { return m_LFpowers; }

  void prepareLFpowers(uint32 startingPoints);

  
 private:
  byte *m_begin;
  uint32 m_length;
  std::vector<uint32> m_LFpowers;
  bool m_isTransformed;

  BWTBlock();
};

} //namespace bwtc

#endif
