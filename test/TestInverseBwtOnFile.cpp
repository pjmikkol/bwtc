/**
 * @file TestInverseBwtOnFile.cpp
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
 * Testing the correctness of BWT inverse transform on file specified by the
 * user.
 */
 
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

int main(int argc, char **argv) {
  int starting_positions = 8;

  if (argc < 2) {
    fprintf(stderr,"usage: %ss filename #starting-positions[default = 8]\n",
        argv[0]);
    exit(1);
  }

  if (argc > 2) {
    int tmp = atoi(argv[2]);
    if (!(tmp >= 1 && tmp <= 256)) {
      fprintf(stderr, "Wrong number of starting positions!\n");
      exit(1);
    }
    starting_positions = tmp;
  }

  FILE *f = fopen(argv[1], "r");
  if (!f) {
    fprintf(stderr,"Cannot open %s\n", argv[1]);
    exit(1);
  }

  fseek(f, 0, SEEK_END);
  long n = ftell(f);
  fseek(f, 0, SEEK_SET);
  byte *t = new byte[n];
  size_t have_read = fread(t, 1, n, f);
  assert(have_read == n);
  fclose(f);
  fprintf(stderr, "File size = %lu bytes.\n", n);
  
  std::vector<byte> v(n + 1);
  std::copy(t, t + n, v.begin());
  BWTransform* transform = giveTransformer('a');
  MainBlock* data;
  std::vector<uint32> LFpowers;
  LFpowers.resize(starting_positions);
  data = new MainBlock(&v, NULL, (uint64)n);
  transform->connect(data);

  fprintf(stderr,"Forward transform... ");
  std::vector<byte>* result = transform->doTransform(LFpowers);
  fprintf(stderr,"DONE\n");

  delete data;
  delete transform;
  delete[] t;
    
  InverseBWTransform* inverse_transform = giveInverseTransformer();
  clock_t time_start = clock();
  fprintf(stderr,"Inverse transform... ");
  std::vector<byte> *original =
    inverse_transform->doTransform(&(*result)[0], result->size(), LFpowers);
  fprintf(stderr,"DONE\n");
  clock_t time_end = clock();

  delete inverse_transform;
  delete result;

  t = new byte[n];
  if (!t) {
    fprintf(stderr, "Cannot allocate memory to hold original file.\n");
    exit(1);
  }
  f = fopen(argv[1], "r");
  if (!f) {
    perror(argv[1]);
    exit(1);
  }
  have_read = fread(t, 1, n, f);
  assert(have_read == n);
  fclose(f);

  bool correct = true;
  if (original->size() != (uint32)n) {
    correct = false;
  }
  for (int j = 0; j < n; ++j) {
    if ((*original)[j] != t[j]) {
      correct = false;
      break;
    }
  }
  if (!correct) {
    fprintf(stderr,"Error.\n");
  } else {
    fprintf(stderr,"Ok.\n");
  }
  fprintf(stderr,"Inversing time: %5.2Lf s\n", ((long double)(time_end - time_start)) / CLOCKS_PER_SEC);

  delete original;
  delete[] t;

  return 0;
}

