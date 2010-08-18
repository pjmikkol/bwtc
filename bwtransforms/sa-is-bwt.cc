#include <cassert>

#include <algorithm>
#include <vector>

#include "sa-is-bwt.h"
#include "../globaldefs.h"


namespace bwtc {

#define CHR(i) ( cs == sizeof(uint32) ? ( (uint32*)s)[i] : ((byte*)s)[i] )
#define isLMS(i) (i > 0 && t[i] && !t[i-1])

/* Finds the start or end of each bucket (depending on end-parameter) */
void SAISBWTransform::GetBuckets(byte *s, uint32 *bkt, uint32 n, uint32 K,
                                 int cs, bool end)
{
  // TODO: fix the handling of K and buckets
  uint32 sum = 0;
  std::fill(bkt, bkt + K + 1, 0);
  for(uint32 i = 0; i < n; ++i) ++bkt[CHR(i)];
  for(uint32 i = 0; i <= K; ++i) {
    sum += bkt[i];
    bkt[i] = end ? sum : sum - bkt[i];
  }
}

void SAISBWTransform::InduceSAl(const std::vector<bool>& t, uint32 *SA, byte *s,
                                uint32 *bkt, uint32 n, uint32 K, int cs,
                                bool end)
{
  GetBuckets(s, bkt, n, K, cs, end);
  for(uint32 i = 0; i < n; ++i) {
    uint32 j = SA[i] - 1;
    if(j < UINT32_MAX - 1 && !t[j]) SA[bkt[CHR(j)]++] = j;
  }
}

void SAISBWTransform::InduceSAs(const std::vector<bool>& t, uint32 *SA, byte *s,
                                uint32 *bkt, uint32 n, uint32 K, int cs,
                                bool end)
{
  GetBuckets(s, bkt, n, K, cs, end);
  for(int64 i = n - 1; i >= 0; --i) {
    uint32 j = SA[i] - 1;
    if(j < UINT32_MAX - 1 && t[j]) SA[--bkt[CHR(j)]] = j;
  }
}

void SAISBWTransform::SA_IS(byte *s, uint32 *SA, uint32 n, uint32 K, int cs)
{
  assert(n >= 2);
  assert(n < UINT32_MAX - 1);

  std::vector<bool> t(n, false);

  /* Classify the characters into L- and S-types. */
  t[n-2] = false;
  t[n-1] = true;
  for(int64 i = n - 3; i >= 0; --i)
    t[i] = ((CHR(i) < CHR(i+1)) || (CHR(i) == CHR(i+1) && t[i+1])) ?
        true : false;

  /* Stage 1: sort all the S-substrings */
  uint32 *bkt = new uint32[K+1];
  GetBuckets(s, bkt, n, K, cs, true);
  std::fill(SA, SA + n, UINT32_MAX);
  for(uint32 i = 1; i < n; ++i)
    if(isLMS(i)) SA[--bkt[CHR(i)]] = i;

  InduceSAl(t, SA, s, bkt, n, K, cs, false);
  InduceSAs(t, SA, s, bkt, n, K, cs, true);

  delete [] bkt;

  /* Compact the sorted substrings into the first n1 items of SA */
  uint32 n1 = 0;
  for(uint32 i = 0; i < n; ++i)
    if(isLMS(SA[i])) SA[n1++] = SA[i];

  assert(2*n1 < n);

  /* Find the lexicographic names of all substrings */
  std::fill(SA + n1, SA + n, UINT32_MAX);
  uint32 name = 0, prev = UINT32_MAX;
  for(uint32 i = 0; i < n1; ++i) {
    uint32 pos = SA[i];
    bool diff = false;
    for(uint32 d = 0; d < n; ++d) {
      if(prev == UINT32_MAX ||CHR(pos+d) != CHR(prev+d)||t[pos+d] != t[prev+d])
      {
        diff = true;
        break;
      } else if (d > 0 && (isLMS(pos+d) || isLMS(prev+d))) {
        break;
      }
    }

    if (diff) {
      ++name;
      prev = pos;
    }
    pos /= 2; //((pos%2 == 0) ? pos:(pos-1)) / 2;
    SA[n1+pos] = name - 1;
  }
  for(uint32 i = n - 1, j = n - 1; i >= n1; --i) {
    if(SA[i] < UINT32_MAX) SA[j--] = SA[i];
  }

  /* Stage 2: Solve the reduced problem. Recurse if names are not unique. */
  uint32 *SA1 = SA, *s1 = SA + n - n1;
  if(name < n1)
    SA_IS((byte*)s1, SA1, n1, name - 1, sizeof(uint32));
  else
    for(uint32 i = 0; i < n1; ++i) SA1[s1[i]] = i;

  /* Stage 3: Induce the result for the original problem. */
  bkt = new uint32[K+1];
  /* Put all LMS-characters into their buckets. */
  GetBuckets(s, bkt, n, K, cs, true);
  for(uint32 i = 1, j = 0; i < n; ++i)
    if(isLMS(i)) s1[j++] = i;
  for(uint32 i = 0; i < n1; ++i) SA1[i] = s1[SA1[i]];
  std::fill(SA + n1, SA + n, UINT32_MAX);
  for(int64 i = n1 - 1; i >= 0; --i) {
    uint32 j = SA[i];
    SA[i] = UINT32_MAX;
    SA[--bkt[CHR(j)]] = j;
  }
  InduceSAl(t, SA, s, bkt, n, K, cs, false);
  InduceSAs(t, SA, s, bkt, n, K, cs, true);

  delete [] bkt;
}

} //namespace bwtc
