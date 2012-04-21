/**
 * @file FrequencyTable.hpp
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
 * Implementation of FrequencyTable-class which is used for bookkeeping
 * the frequencies of bytes.
 */

#ifndef FREQUENCY_TABLE_HPP_
#define FREQUENCY_TABLE_HPP_

#include "../globaldefs.hpp"

#include <utility>

namespace bwtc {

/* Data structure for holding the frequencies of bytes. */
class FrequencyTable {
 public:
  FrequencyTable();
  FrequencyTable(const FrequencyTable& freqTable);

  FrequencyTable& operator=(const FrequencyTable& freqTable);

  /** Constructs FrequencyTable from given freqs */
  FrequencyTable(uint32* frequencies); 

  /**Does basically the same as the above ctor. */
  void initialize(uint32* frequencies);

  /** Returns the i:th lowest freq*/
  uint32 getFrequency(int i) const; 

  uint32 getFrequencyWithKey(byte key) const;

  /** Returns the key which has i:th lowest freq*/
  byte getKey(int i) const; 
  
  /**Decrease value of given key.*/
  bool decrease(byte key, uint32 decrement);
  
  /**Increase value of given key.*/
  void increase(byte key, uint32 increment);

 private:
  void initLocations();
  bool test();
  std::pair<byte, uint32> m_frequencies[256];
  byte m_location[256];
};


} //namespace bwtc

#endif
