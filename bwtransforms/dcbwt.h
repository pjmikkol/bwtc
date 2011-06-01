/**************************************************************************
 *  Copyright 2007, Google Inc.                                           *
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

/***********************************************************************
 * Burrows-Wheeler transform using the Karkkainen-Burkhardt algorithm  *
 * which is based on difference cover sampling.                        *
 *                                                                     *
 * Original implementation of the algorithm can be found at:           *
 * http://code.google.com/p/dcs-bwt-compressor/                        *
 ***********************************************************************/
#ifndef DCBW_TRANSFORM_H_
#define DCBW_TRANSFORM_H_

#include "../MainBlock.hpp"
#include "bw_transform.h"

namespace bwtc {

class DCBWTransform : public BWTransform {
 public:
  DCBWTransform(int period);
  virtual ~DCBWTransform();
  
  virtual std::vector<byte>* DoTransform(uint64* eob_byte);

  virtual uint64 MaxSizeInBytes(uint64 block_size) const;
  virtual uint64 MaxBlockSize(uint64 memory_budget) const;
  virtual uint64 SuggestedBlockSize(uint64 memory_budget) const;

 private:
  static const uint64 kMemoryOverhead = (1 << 20);
  static const int kMinLogPeriod = 3;
  static const int kMaxLogPeriod = 15;
  int log_period_;
  uint32 period_;
};

} //namespace bwtc

#endif
