
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

#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

#include "../MainBlock.hpp"
#include "../bwtransforms/bw_transform.h"
#include "../bwtransforms/dcbwt.h"

namespace bwtc {
int verbosity = 7;
}

using bwtc::uint64;
using bwtc::byte;

namespace tests {

void SimpleBWTtest(const char* input, uint64 length) {
  bwtc::BWTransform* bwt = bwtc::GiveTransformer('d');
  std::vector<byte> data(&input[0], &input[length]);
  std::vector<uint64> dummy;
  bwtc::MainBlock block(&data, &dummy, length);
  bwt->Connect(&block);
  std::reverse(data.begin(), data.end());
  uint64 eob;
  std::vector<byte>* result = bwt->DoTransform(&eob);
  for(unsigned i = 0; i < length + 1; ++i) {
    if(i == eob) std::cout << '$';
    else std::cout << (*result)[i];
  }
  std::cout << "\n" << "eob-byte position : " << eob << "\n";
  delete bwt;
  delete result;
}

void BWTfromFile(const char *fname) {
  bwtc::BWTransform* bwt = bwtc::GiveTransformer('d');
  std::ifstream file(fname);
  std::vector<byte> data((std::istream_iterator<byte>(file)),
                         std::istream_iterator<byte>());
  std::vector<uint64> dummy;
  uint64 length = data.size();
  bwtc::MainBlock block(&data, &dummy, length);
  bwt->Connect(&block);
  std::reverse(data.begin(), data.end());
  uint64 eob;
  std::vector<byte>* result = bwt->DoTransform(&eob);
  for(unsigned i = 0; i < length + 1; ++i) {
    if(i == eob) std::cout << '$';
    else std::cout << (*result)[i];
  }
  delete bwt;
  delete result;
}


} //namespace tests

int main(int argc, char **argv) {
  std::string tst("mississippi");
  tests::SimpleBWTtest(tst.c_str(), tst.size());
  if (argc > 1) tests::BWTfromFile(argv[1]);
}
