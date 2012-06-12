/**
 * @file Decompressor.cpp
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
 * Implementation of Decompressor-class.
 */

#include "Decompressor.hpp"
#include "Streams.hpp"
#include "preprocessors/Postprocessor.hpp"

#include <string>

namespace bwtc {

Decompressor::Decompressor(const std::string& in, const std::string& out)
    : m_in(new InStream(in)), m_out(new OutStream(out)),
      m_decoder(0), m_postprocessor(verbosity > 1) {}

Decompressor::Decompressor(InStream* in, OutStream* out)
    : m_in(in), m_out(out), m_decoder(0), m_postprocessor(verbosity > 1) {}

Decompressor::~Decompressor() {
  delete m_in;
  delete m_out;
  delete m_decoder;
}

size_t Decompressor::readGlobalHeader() {
  char entropyDecoder = static_cast<char>(m_in->readByte());
  delete m_decoder;
  m_decoder = giveEntropyDecoder(entropyDecoder);
  return 1;
}

} //namespace bwtc
