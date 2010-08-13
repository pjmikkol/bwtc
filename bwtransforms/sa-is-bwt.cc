#include <algorithm>
#include <vector>

#include "sa-is-bwt.h"

namespace bwtc {

#define chr(i) ( cs == sizeof(int) ? ( (int*)s)[i] : ((byte*)s)[i] )
#define isLMS(i) (i > 0 && t[i] && !t[i-1])

void SAISBWTransform::GetBuckets(
    byte *s, int *bkt, int n, int K, int cs, bool end)
{
  // TODO: fix the handling of K and buckets
  int sum = 0;
  std::fill(bkt, bkt + K + 1, 0);
  for(int i = 0; i < n; ++i) ++bkt[chr(i)];
  for(int i = 0; i <= K; ++i) {
    sum += bkt[i];
    bkt[i] = end ? sum : sum - bkt[i];
  }
}

void SAISBWTransform::InduceSAl(std::vector<bool> *t, int *SA, byte *s,
                                int *bkt, int n, int K, int cs, bool end)
{
  GetBuckets(s, bkt, n, K, cs, end);
  for(int i = 0; i < n; ++i) {
    int j = SA[i] - 1;
    if(j >= 0 && !(*t)[j]) SA[bkt[chr(j)]++] = j;
  }
}

void SAISBWTransform::InduceSAs(std::vector<bool> *t, int *SA, byte *s,
                                int *bkt, int n, int K, int cs, bool end)
{
  GetBuckets(s, bkt, n, K, cs, end);
  for(int i = n - 1; i >= 0; --i) {
    int j = SA[i] - 1;
    if(j >= 0 && (*t)[j]) SA[--bkt[chr(j)]] = j;
  }
}

void SAISBWTransform::SA_IS(byte *s, int *SA, int n, int K, int cs)
{
  std::vector<bool> t(n, false);

  t[n-2] = false;
  t[n-1] = true;
  for(int i = n - 3; i >= 0; --i)
    t[i] = ((s[i] < s[i+1]) || (s[i] == s[i+1] && t[i+1])) ? true : false;

  int *bkt = new int[K+1];
  GetBuckets(s, bkt, n, K, cs, true);
  std::fill(SA, SA + n, -1);
  //TODO: What happens to sentinel?
  for(int i = 1; i < n; ++i)
    if(isLMS(i)) SA[--bkt[chr(i)]] = i;

  InduceSAl(&t, SA, s, bkt, n, K, cs, false);
  InduceSAs(&t, SA, s, bkt, n, K, cs, true);

  delete [] bkt;

  int n1 = 0;
  for(int i = 0; i < n; ++i)
    if(isLMS(SA[i])) SA[n1++] = SA[i];

  std::fill(SA + n1, SA + n, -1);
  int name = 0, prev = -1;
  for(int i = 0; i < n1; ++i) {
    int pos = SA[i];
    bool diff = false;
    for(int d = 0; d < n; ++d) {
      if(prev == -1||chr(pos+d) != chr(prev+d)||t[pos+d] != t[prev+d]) {
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
    pos = ((pos%2 == 0) ? pos:(pos-1)) / 2;
    SA[n1+pos] = name - 1;
  }

  for(int i = n - 1, j = n - 1; i >= n1; --i) {
    if(SA[i] >= 0) SA[j--] = SA[i];
  }

  int *SA1 = SA, *s1 = SA + n - n1;
  if(name < n1)
    SA_IS((byte*)s1, SA1, n1, name - 1, sizeof(int));
  else
    for(int i = 0; i < n1; ++i) SA1[s1[i]] = i;

  bkt = new int[K+1];
  GetBuckets(s, bkt, n, K, cs, true);
  for(int i = 1, j = 0; i < n; ++i)
    if(isLMS(i)) s1[j++] = i;
  for(int i = 0; i < n1; ++i) SA1[i] = s1[SA1[i]];
  std::fill(SA + n1, SA + n, -1);
  for(int i = n1 - 1; i >= 0; --i) {
    int j = SA[i];
    SA[i] = -1;
    SA[--bkt[chr(j)]] = j;
  }
  InduceSAl(&t, SA, s, bkt, n, K, cs, false);
  InduceSAs(&t, SA, s, bkt, n, K, cs, true);

  delete [] bkt;
}

} //namespace bwtc
