/**
 * @file DivsufsortTest.cpp
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
 * Test of the libivsufsort.
 */

#include <cstdlib>
#include <cstring>
#include <ctime>

#include <iostream>
#include <vector>

#include "../globaldefs.hpp"
#include "../bwtransforms/divsufsort.h"
#include "../MainBlock.hpp"

#undef NDEBUG

//using bwtc::uint64;
using bwtc::uint32;
using bwtc::byte;

namespace tests {

void BwtSaisTest(char *arg) {
  int len = strlen(arg);
  byte *str = new byte[len+1];
  byte *res = new byte[len+1];
  strcpy((char*)str, arg);

  std::vector<uint32> LFpowers;
  LFpowers.resize(1);
  divbwt(str, res, 0, len+1, &LFpowers[0], LFpowers.size());

  for(int i = 0; i < len + 1; ++i) {
    //std::cout << res[i];
    if(i != LFpowers[0]) std::cout << res[i];
    else std::cout << '$';
  }
  std::cout << "\nReturn value: " << LFpowers[0] << "\n";
  delete [] str;
  delete [] res;
}

}


int main(int argc, char** argv) {
  using namespace tests;
  if(argc > 1) {
    BwtSaisTest(argv[1]);
  }
}
