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

#ifndef BWTC_TEST_PREPROCESSOR_H_
#define BWTC_TEST_PREPROCESSOR_H_

#include <iostream> /* for std::streamsize*/
#include <string>

#include "../block.h"
#include "../block_manager.h"
#include "../globaldefs.h"
#include "../stream.h"
#include "preprocessor.h"

namespace bwtc {
/* This class is meant for testing the different preprocessor-options.
 * It can be done by using the available public-functions such as
 * CompressPairs */
class TestPreProcessor : public PreProcessor {
 public:
  TestPreProcessor(uint64 block_size);
  virtual ~TestPreProcessor();
  /* Reads and preprocesses data to byte array provided by block_manager_*/
  uint64 CompressPairs();
  /* Same as above, but different compression algorithm */
  uint64 CompressRuns();
  /* Initialize target-array for reading */
  void InitializeTarget();
  /* Fills the buffer from instream. returns true if something is read*/
  uint64 FillBuffer();

  MainBlock *curr_block_;
};

} // namespace bwtc

#endif
