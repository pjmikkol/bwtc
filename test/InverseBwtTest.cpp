/**
 * @file InverseBwtTest.cpp
 * @author Dominik Kempa <dominik.kempa@cs.helsinki.fi>
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
 * Testing the correctness of BWT inverse transform on many random inputs.
 */

#define MAIN

#include <cstdio>
#include <cstring>
#include <string>
#include <ctime>
#include <cstdlib>
#include <unistd.h>
#include <algorithm>
#include <cassert>
#include <vector>

#include "../globaldefs.hpp"
#include "../bwtransforms/BWTransform.hpp"
#include "../bwtransforms/InverseBWT.hpp"
#include "../bwtransforms/MtlSaInverseBWT.hpp"
#include "../bwtransforms/SA-IS-bwt.hpp"

using namespace bwtc;

inline int my_random(int p, int r) {
  return p + rand() % (r - p + 1);
}

void test(byte *t, uint32 n) {
  std::vector<byte> v(n + 1);
  std::copy(t, t + n, v.begin());
  BWTransform* transform = giveTransformer('a');
  InverseBWTransform* inverse_transform = giveInverseTransformer();

  std::vector<byte> data(t, t+n);
  data.push_back(0);

  std::vector<uint32> LFpowers;
  int starting_points = my_random(1, std::min(256, (int)n));
  LFpowers.resize(starting_points);

  //std::vector<byte>* result =
  transform->doTransform(&data[0], n+1, LFpowers);

  std::vector<byte> *original =
      inverse_transform->doTransform(&data[0], n, LFpowers);

  bool ok = true;
  if (original->size() != n) ok = false;
  for (uint32 j = 0; j < n; ++j) {
    if ((*original)[j] != t[j]) {
      ok = false;
      break;
    }
  }
  if (!ok) {
    fprintf(stderr,"FAIL\n");
    fprintf(stderr,"t = ");
    for (uint32 j = 0; j < n; ++j) {
      fprintf(stderr,"%c", t[j]);
    }
    fprintf(stderr,"\n");
    fprintf(stderr,"returned = ");
    for (uint32 j = 0; j < n; ++j) {
      fprintf(stderr,"%c", data[j]);
    }
    fprintf(stderr,"\n");
    exit(1);
  }

  delete original;
  delete transform;
  delete inverse_transform;
}

void test_random(int testcases, int max_n, int max_sigma) {  
  fprintf(stderr,"TEST, testcases = %d, max_n = %d, max_sigma = %d\n",
      testcases, max_n, max_sigma);
  int progress_step = testcases / 1000;
  if (progress_step == 0) {
    progress_step = 1;
  }
  for (int tc = 0; tc < testcases; ++tc) {
    if ( ( tc % progress_step ) == 0 ) {
      fprintf(stderr,"Progress: %5.2Lf%%\r",
          ((long double) 100.0 * tc ) / testcases);
    }
    uint32 n = my_random(1, max_n);
    byte *t = new byte[n];
    for (uint32 j = 0; j < n; ++j) {
      t[j] = my_random(0, max_sigma - 1);
    }
    test(t, n);
    delete[] t;
  }
}

int main() {
  srand(time(0) + getpid());
  test_random(1000000,    3, 256);
  test_random(1000000,   10, 256);
  test_random(100000,   100, 256);
  test_random(100000,  1000, 256);
  test_random(10000,  10000, 256);
  test_random(10000, 100000, 256);
  test_random(1000, 1000000, 256);
  fprintf(stderr,"All tests passed.\n");
  return 0;
}

