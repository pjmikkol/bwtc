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
  MainBlock* data;
  std::vector<uint32> LFpowers;
  int starting_points = my_random(1, 256);
  LFpowers.resize(starting_points);
  data = new MainBlock(&v, NULL, (uint64)n);
  transform->connect(data);
  std::vector<byte>* result = transform->doTransform(LFpowers);
  std::vector<byte> *original =
    inverse_transform->doTransform(&(*result)[0], result->size(), LFpowers);

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
      fprintf(stderr,"%c", (*result)[j]);
    }
    fprintf(stderr,"\n");
    exit(1);
  }

  delete original;
  delete transform;
  delete inverse_transform;
  delete result;
  delete data;
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

