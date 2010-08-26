/**************************************************************************
 *  Copyright 2010, Pekka Mikkola, pjmikkol (at) cs.helsinki.fi           *
 *                                                                        *
 *  This file is part of bwtc.                                            *
 *                                                                        *
 *  bwtc is free software: you can redistribute it and/or modify          *
 *  it under the terms of the GNU General Public License as published by  *
 *  the Free Software Foundation, either version 3 of the License, or     *
 *  (at your option) any later version.                                   *
 *                                                                        *
 *  bwtc is distributed in the hope that it will be useful,               *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *  GNU General Public License for more details.                          *
 *                                                                        *
 *  You should have received a copy of the GNU General Public License     *
 *  along with bwtc.  If not, see <http://www.gnu.org/licenses/>.         *
 **************************************************************************/

#ifndef BWTC_PREPROCESSOR_H_
#define BWTC_PREPROCESSOR_H_

#include <iostream> /* for std::streamsize*/
#include <string>

#include "../block.h"
#include "../block_manager.h"
#include "../globaldefs.h"
#include "../stream.h"

namespace bwtc {

/* Floor of logarithm of base two */
byte LogFloor(unsigned n);

unsigned MostSignificantBit16(unsigned n);

unsigned MostSignificantBit(unsigned n);

//TODO: Make this class to abstract class (interface for real implementations)
class PreProcessor {
 public:
  PreProcessor(uint64 block_size);
  virtual ~PreProcessor();
  virtual void Connect(std::string source_name);
  virtual void AddBlockManager(BlockManager* bm);
  /* Reads and preprocesses data to byte array provided by block_manager_*/
  virtual MainBlock* ReadBlock();

 protected:
  InStream* source_;
  uint64 block_size_;
  BlockManager* block_manager_;

 private:
  /* This should be done during preprocessing*/
  void BuildStats(std::vector<byte>* data, std::vector<uint64>* stats,
                  uint64 data_size);
  PreProcessor& operator=(const PreProcessor& p);
  PreProcessor(const PreProcessor&);
};

/* Data structure for holding the frequencies of bytes. */
class FreqTable {
 public:
  FreqTable();
  FreqTable(uint64* frequencies); /* Constructs FreqTable from given freqs */
  const uint64& operator[](unsigned i); /* Returns the i:th lowest freq */
  byte Key(unsigned i); /* Returns the key which has i:th lowest freq */
  bool Decrease(unsigned key, uint64 decrement);
  void Increase(unsigned key, uint64 increment);

 private:
  void InitLocations();
  bool Test();
  std::pair<byte, uint64> freq_[256];
  byte location_[256];
};

/* This function returns chosen preprocessor */ 
PreProcessor* GivePreProcessor(
    char choice, uint64 block_size, const std::string& input);

uint64 CompressCommonPairs(byte *from, uint64 length);
uint64 CompressLongRuns(byte *from, uint64 length);

} // namespace bwtc


#endif
