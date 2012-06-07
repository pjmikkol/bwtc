/**
 * @file sa-is_test.cc
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
 * Test of the SA-IS BWT-algorithm implementations.
 */

#include <cstdlib>
#include <cstring>
#include <ctime>

#include <iostream>

#include "../globaldefs.hpp"
#include "../bwtransforms/SA-IS-bwt.hpp"
#include "../bwtransforms/sais.hxx"

#undef NDEBUG

//using bwtc::uint64;
using bwtc::uint32;
using bwtc::byte;

namespace tests {

int SufCmp(byte *str, uint32 s1, uint32 s2, unsigned n) {
  unsigned limit = n - ((s1 > s2)? s1 : s2);
  for(unsigned i = 0; i < limit; ++i) {
    if(str[s1 + i] < str[s2 + i]) return -1;
    else if (str[s1 + i] > str[s2 + i]) return 1;
  }
  if(s1 > s2) return -1;
  else return 1;
}

void test_sais() {
  int size = (rand() & 0x0000FFFF) + 1;
  byte *str = new byte[size+1];
  for(int i = 0; i < size; ++i) {
    str[i] = rand() & 0xFF;
  }
  str[size] = 0;
  int *SA = new int[size + 1];
  saisxx(str, SA, size + 1, 256);
  for(int i = 1; i < size; ++i) {
    assert(SufCmp(str, SA[i], SA[i+1], size) < 0);
    assert(SA[i] < size);
  }
  std::cout << "." << std::flush;
  assert(SA[0] == size);
}


void BwtSaisTest(char *arg) {
  int len = strlen(arg);
  byte *str = new byte[len+1];
  strcpy((char*)str, arg);
  int *SA = new int[len+1];
  byte *res = new byte[len+1];
  std::vector<uint32> LFpowers;
  LFpowers.resize(1);
  saisxx_bwt(str, res, SA, len + 1, LFpowers, 256);
  int val = LFpowers[0];
  for(int i = 0; i < len + 1; ++i) {
    if(i != val) std::cout << res[i];
    else std::cout << '$';
  }
  std::cout << "\nReturn value: " << val << "\n";
  delete [] SA;
  delete [] str;
}

}


int main(int argc, char** argv) {
  using namespace tests;
  if(argc > 1) {
    BwtSaisTest(argv[1]);
  }
  srand( time(NULL));
  for(int i = 0; i < 10; ++i) {
    test_sais();
  }
  std::cout << "\nSA_IS passed all tests.\n";
}
