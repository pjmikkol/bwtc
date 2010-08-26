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

#ifndef BWTC_POSTPROCESSOR_H_
#define BWTC_POSTPROCESSOR_H_

#include <vector>

#include "../globaldefs.h"

//TODO: PostProcessor-class since aAt the moment only algorithms for
//      postprocessing exists

namespace bwtc {
  
uint64 UncompressCommonPairs(std::vector<byte> *from, uint64 length);
uint64 UncompressLongRuns(std::vector<byte> *from, uint64 length);
uint64 UncompressSequences(std::vector<byte> *from, uint64 length);

}

#endif
