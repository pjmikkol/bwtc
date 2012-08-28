/**
 * @file BWTManager.cpp
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
 * Implementation of BWT-manager.
 */

#include "../globaldefs.hpp"
#include "BWTManager.hpp"
#include "BWTransform.hpp"
#include "SA-IS-bwt.hpp"
#include "Divsufsorter.hpp"

namespace bwtc {

BWTManager::BWTManager() : m_startingPoints(1) {}

BWTManager::BWTManager(uint32 startingPoints)
    : m_startingPoints(startingPoints) {}

BWTManager::~BWTManager() {
  for(size_t i = 0; i < m_transformers.size(); ++i) {
    delete m_transformers[i];
  }
}

void BWTManager::doTransform(BWTBlock& block) {
  assert(!block.isTransformed());
  block.prepareLFpowers(m_startingPoints);
  //Something more sophisticated here if choosing algorithm automatically:
  m_transformers[0]->doTransform(block);
}

void BWTManager::doTransform(BWTBlock& block, uint32 *freqs) {
  assert(!block.isTransformed());
  block.prepareLFpowers(m_startingPoints);
  //Something more sophisticated here if choosing algorithm automatically:
  m_transformers[0]->doTransform(block, freqs);
}

void BWTManager::setStartingPoints(uint32 startingPoints) {
  if(startingPoints < 1) startingPoints = 1;
  else if(startingPoints > 256) startingPoints = 256;
  m_startingPoints = startingPoints;
}

uint32 BWTManager::getStartingPoints() const {
  return m_startingPoints;
}

bool BWTManager::isValidChoice(char c) {
  return c == 'd' || c == 's' || c == 'a';
}

void BWTManager::initialize(char choice) {
  if(choice == 's') {
    m_transformers.push_back(new SAISBWTransform());
  } else {
    m_transformers.push_back(new Divsufsorter());
  }
}

} //namespace bwtc
