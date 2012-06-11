/**
 * @file Precompressor.hpp
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
 * Header for precompressor.
 */

#ifndef BWTC_PRECOMPRESSOR_HPP_
#define BWTC_PRECOMPRESSOR_HPP_

#include <iostream> /* for std::streamsize*/
#include <string>

#include "../globaldefs.hpp"
#include "../Streams.hpp"
#include "../PrecompressorBlock.hpp"

namespace bwtc {

class Precompressor {
 public:
  Precompressor();
  Precompressor(const std::string& prepr);
  ~Precompressor();

  /* Reads and preprocesses data to byte array */
  PrecompressorBlock* readBlock(size_t blockSize, RawInStream* in) const;

  // TODO: give also the temporary memory to parameter (to be able
  // precompress in-place)
  void precompress(PrecompressorBlock& block) const;

  const std::string options() const { return m_preprocessingOptions; }
  
 private:
  Precompressor& operator=(const Precompressor& p);
  Precompressor(const Precompressor&);

  const std::string m_preprocessingOptions;
};

} // namespace bwtc


#endif
