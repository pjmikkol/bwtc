/**
 * @file Compressor.cpp
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
 * Implementation of Compressor.
 */

#include "Compressor.hpp"

#include <string>

namespace bwtc {

Compressor::
Compressor(const std::string& in, const std::string& out, size_t memLimit)
    : m_in(new RawInStream(in)), m_out(new RawOutStream(out)), m_coder(0),
      m_preprocessor(0), m_memLimit(memLimit) {}

Compressor::Compressor(RawInStream* in, RawOutStream* out, size_t memLimit)
    : m_in(in), m_out(out), m_coder(0), m_preprocessor(0),
      m_memLimit(memLimit) {}

Compressor::~Compressor() {
  delete m_in;
  delete m_out;
  delete m_coder;
  delete m_preprocessor;
}

void Compressor::setEntropyEncoder(char choice) {
  delete m_coder;
  //m_coder = giveEntropyEncoder(choice);
}

} //namespace bwtc
