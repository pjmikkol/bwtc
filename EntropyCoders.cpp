/**
 * @file EntropCoders.cpp
 * @author Dominik Kempa <dominik.kempa@cs.helsinki.fi>
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
 * Implementation of base classes for entropy coders.
 *
 */

#include <iostream>
#include <string>
#include <vector>

#include "EntropyCoders.hpp"
#include "Coders.hpp"
#include "WaveletCoders.hpp"
#include "HuffmanCoders.hpp"

namespace bwtc {

EntropyEncoder* giveEntropyEncoder(const std::string& destination,
    char prob_model) {
  return new HuffmanEncoder(destination, prob_model);
  //return new WaveletEncoder(destination, prob_model);
}

EntropyDecoder* giveEntropyDecoder(const std::string& source) {
 return new HuffmanDecoder(source);
 //return new WaveletDecoder(source);
}

} // namespace bwtc
