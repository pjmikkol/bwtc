/**
 * @file Postprocessor.hpp
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
 * Header for the inverses of preprocessing algorithms.
 */

#ifndef BWTC_POSTPROCESSOR_HPP_
#define BWTC_POSTPROCESSOR_HPP_

#include <vector>

#include "../globaldefs.hpp"

//TODO: PostProcessor-class since at the moment only algorithms for
//      postprocessing exists

namespace bwtc {

class PostProcessor {
 public:
  struct Replacement {
    Replacement();
    Replacement(uint32 len, uint16 repl, bool pair);
    Replacement(const Replacement& r);
    Replacement& operator=(const Replacement& repl);
    
    uint32 length;
    uint16 replacement;
    bool isPair;
  };

  PostProcessor(const std::string& postProcOptions);

  void postProcess(std::vector<byte> *data);

  static size_t uncompressCommonPairs(std::vector<byte> *from, size_t length);
  static size_t uncompressLongRuns(std::vector<byte> *from, size_t length);
  static size_t uncompressSequences(std::vector<byte> *from, size_t length);
  static size_t uncompressPairsAndRuns(std::vector<byte> *compressed, size_t length);

 private:
  std::string m_options;
};
  
} //namespace bwtc

#endif
