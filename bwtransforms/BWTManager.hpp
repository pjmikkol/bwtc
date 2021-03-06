/**
 * @file BWTManager.hpp
 * @author Pekka Mikkola <pmikkol@gmail.com>
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
 * Header for BWT-manager. The choice of BWT-algorithm is done in this
 * class. In future, if needed, the choice of the algorithm can be done in
 * fly based on the running times.
 *
 */

#ifndef BWTC_BWTMANAGER_HPP_
#define BWTC_BWTMANAGER_HPP_

#include "../globaldefs.hpp"
#include "../BWTBlock.hpp"
#include "BWTransform.hpp"

#include <vector>

namespace bwtc {

class BWTManager {
 public:
  BWTManager();
  BWTManager(uint32 startingPoints);
  ~BWTManager();

  void doTransform(BWTBlock& block);
  void doTransform(BWTBlock& block, uint32 *freqs);
  void initialize(char choice);
  void setStartingPoints(uint32 startingPoints);
  uint32 getStartingPoints() const;

  static bool isValidChoice(char c);
  
 private:
  std::vector<BWTransform*> m_transformers;
  uint32 m_startingPoints;
};

}  //namespace bwtc

#endif
