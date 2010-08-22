/*******************************************************************
 * Burrows-Wheeler transform using the SA-IS algorithm desribed in *
 * "Two Efficient Algorithms for Linear Suffix Array Construction" *
 * by Nong, Zhang & Chan                                           *
 ******************************************************************/
#ifndef BWTC_SAIS_BWT_H_
#define BWTC_SAIS_BWT_H_

#include <cassert>

#include <algorithm>
#include <vector>



#include <iostream>




#include "bw_transform.h"
#include "../globaldefs.h"

namespace bwtc {

class SAISBWTransform : public BWTransform {
 public:
  SAISBWTransform();
  virtual ~SAISBWTransform() {}
  virtual std::vector<byte>* DoTransform(uint64* eob_byte);

  /* The following values aren't correct */
  virtual uint64 MaxSizeInBytes(uint64 block_size) const { return 0; }
  virtual uint64 MaxBlockSize(uint64 memory_budget) const { return 0; }
  virtual uint64 SuggestedBlockSize(uint64 memory_budget) const { return 0; }
  
};

} //namespace bwtc


namespace sa_is {

/****************************************************************************
 * We need two slightly different versions of main algorithm. If zero can   *
 * be found from the given string then one have to use SA_IS_Zero.. function*
 * T is the type of the given string and S is the type of the suffix array. *
 * If string is long (4294967295 or longer) then S needs to be uint64.      *
 * Arguments: s is string, SA is suffix array, n is the length of SA and s, *
 * and K is the size of alphabet.                                           *
 * Both versions require that last character of the input string is 0.      *
 ****************************************************************************/
template <typename T, typename S>
void SA_IS_ZeroInclude(T *s, S *SA, int64 n, uint32 K);
template <typename T, typename S>
void SA_IS(T *s, S *SA, int64 n, uint32 K);

namespace {
/* Empty namespace is for the helper procedures. Which are not meant to be *
 * used by themselves. */
template <typename T>
void GetBuckets(T *s, uint32 *bkt, int64 n, uint32 K, bool end);

template <typename T>
void GetBucketsZeroInclude(T *s, uint32 *bkt, int64 n, uint32 K, bool end);

/* Mutual part for the two SA_IS-procedures */
template <typename T, typename S>
void MutualPart(T* s, S *SA, int64 n, uint32 K, std::vector<bool> *isS);

template <typename T>
void ClassifyCharacters(T *s, std::vector<bool>* t, int64 len);

template <typename T, typename S>
void InduceSAl(const std::vector<bool>& t, S *SA, T *s, uint32 *bkt,
               int64 n, uint32 K, bool end);

template <typename T, typename S>
void InduceSAs(const std::vector<bool>& t, S *SA, T *s, uint32 *bkt,
               int64 n, uint32 K, bool end);

} //namespace

} //namespace sa_is


/* Implementation for sa_is-algorithm */
namespace sa_is {

#define isLMS(i) (i < Max<S>::max && i > 0 && t[i] && !t[i-1])

namespace {

template <typename T, typename S>
void MutualPart(T* s, S *SA, int64 n, uint32 K, std::vector<bool> *isS)
{
  std::vector<bool>& t = *isS;
  /* Compact the sorted substrings into the first n1 items of SA */
  int64 n1 = 0;
  for(int64 i = 0; i < n; ++i) 
    if(isLMS(SA[i])) SA[n1++] = SA[i];

  assert(2*n1 <= n);

  /* Find the lexicographic names of all substrings */
  std::fill(SA + n1, SA + n, Max<S>::max);
  S name = 0, prev = Max<S>::max;
  for(int64 i = 0; i < n1; ++i) {
    S pos = SA[i];
    bool diff = false;
    for(S d = 0; d < n; ++d) {
      if(prev == Max<S>::max ||s[pos+d] != s[prev+d]||t[pos+d] != t[prev+d])
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
    pos /= 2; 
    SA[n1+pos] = name - 1;
  }
  for(int64 i = n - 1, j = n - 1; i >= n1; --i) {
    if(SA[i] < Max<S>::max) SA[j--] = SA[i];
  }

  /* Stage 2: Solve the reduced problem. Recurse if names are not unique. */
  S *SA1 = SA, *s1 = SA + n - n1;
  if(name < n1)
    SA_IS(s1, SA1, n1, name - 1);
  else
    for(int64 i = 0; i < n1; ++i) SA1[s1[i]] = i;

  /* Stage 3: Induce the result for the original problem. */
  uint32 *bkt = new uint32[K+1];
  /* Put all LMS-characters into their buckets. */
  GetBuckets(s, bkt, n, K, true);
  for(int64 i = 1, j = 0; i < n; ++i) 
    if(isLMS(i)) s1[j++] = i; 
  for(int64 i = 0; i < n1; ++i) SA1[i] = s1[SA1[i]];
  std::fill(SA + n1, SA + n, Max<S>::max);
  for(int64 i = n1 - 1; i >= 0; --i) {
    S j = SA[i];
    SA[i] = Max<S>::max;
    SA[--bkt[s[j]]] = j;
  }
  InduceSAl(t, SA, s, bkt, n, K, false);
  InduceSAs(t, SA, s, bkt, n, K, true);

  delete [] bkt;
}

/* Finds the start or end of each bucket (depending on end-parameter) */
template <typename T>
void GetBuckets(T *s, uint32 *bkt, int64 n, uint32 K, bool end)
{
  int64 sum = 0;
  std::fill(bkt, bkt + K + 1, 0);
  for(int64 i = 0; i < n; ++i) ++bkt[s[i]];
  for(int64 i = 0; i <= K; ++i) {
    sum += bkt[i];
    bkt[i] = end ? sum : sum - bkt[i];
  }
}

template <typename T>
void GetBucketsZeroInclude(T *s, uint32 *bkt, int64 n, uint32 K, bool end)
{
  int64 sum = 0;
  std::fill(bkt, bkt + K + 1, 0);
  for(int64 i = 0; i < n - 1; ++i) ++bkt[s[i] + 1];
  bkt[0] = 1;
  for(int64 i = 0; i <= K; ++i) {
    sum += bkt[i];
    bkt[i] = end ? sum : sum - bkt[i];
  }
}

template <typename T>
void ClassifyCharacters(T *s, std::vector<bool> *isS, int64 n) {
  std::vector<bool>& t = *isS;
  t[n-2] = false;
  t[n-1] = true;
  for(int64 i = static_cast<int64>(n) - 3; i >= 0; --i)
    t[i] = ((s[i] < s[i+1]) || (s[i] == s[i+1] && t[i+1])) ?
        true : false;
}

template <typename T, typename S>
void InduceSAl(const std::vector<bool>& t, S *SA, T *s, uint32 *bkt, int64 n,
               uint32 K, bool end)
{
  GetBuckets(s, bkt, n, K, end);
  for(int64 i = 0; i < n; ++i) {
    S j = SA[i] - 1;
    if(j < Max<S>::max - 1 && !t[j]) SA[bkt[s[j]]++] = j;
  }
}

template <typename T, typename S>
void InduceSAs(const std::vector<bool>& t, S *SA, T *s, uint32 *bkt, int64 n,
               uint32 K, bool end)
{
  GetBuckets(s, bkt, n, K, end);
  for(int64 i = n - 1; i >= 0; --i) {
    S j = SA[i] - 1;
    if(j < Max<S>::max - 1 && t[j]) SA[--bkt[s[j]]] = j;
  }
}

} //namespace

template <typename T, typename S>
void SA_IS(T *s, S *SA, int64 n, uint32 K)
{
  assert(n >= 2);
  assert(n < Max<S>::max - 1);
  //assert(s[n-1] == 0);

  std::vector<bool> t(n, false);

  /* Classify the characters into L- and S-types. */
  ClassifyCharacters(s, &t, n);

  /* Stage 1: sort all the S-substrings */
  uint32 *bkt = new uint32[K+1];
  GetBuckets(s, bkt, n, K, true);
  std::fill(SA, SA + n, Max<S>::max);
  for(int64 i = 1; i < n; ++i)
    if(isLMS(i)) SA[--bkt[s[i]]] = i;

  InduceSAl(t, SA, s, bkt, n, K, false);
  InduceSAs(t, SA, s, bkt, n, K, true);

  delete [] bkt;
  MutualPart(s, SA, n, K, &t);
}

template <typename T, typename S>
void SA_IS_ZeroInclude(T *s, S *SA, int64 n, uint32 K)
{
  assert(n >= 2);
  assert(n < Max<S>::max - 1);
  assert(s[n-1] == 0);

  std::vector<bool> t(n, false);
  /* Classify the characters into L- and S-types. */
  ClassifyCharacters(s, &t, n);

  /* Stage 1: sort all the S-substrings */
  uint32 *bkt = new uint32[K+1];
  GetBucketsZeroInclude(s, bkt, n, K, true);
  std::fill(SA, SA + n, Max<S>::max);
  for(int64 i = 1; i < n - 1; ++i)
    if(isLMS(i)) SA[--bkt[s[i] + 1]] = i;
  SA[--bkt[0]] = n - 1;
  InduceSAl(t, SA, s, bkt, n, K, false);
  InduceSAs(t, SA, s, bkt, n, K, true);

  delete [] bkt;
  MutualPart(s, SA, n, K, &t);
  SA[0] = n - 1;
}

namespace {
const uint32 kMinBucketSize = 256;
const uint32 kSaisLMSsort2Limit = 0x3fffffff;

template <typename T, typename S>
void getCounts(const T *s, S *C, S n, S K) {
  std::fill(C, C + K, 0);
  for(S i = 0; i < n; ++i) ++C[s[i]];
}

template <typename S>
void getBuckets(const S *C, S *B, S K, bool end) {
  S sum = 0;
  if(end) {
    for(S i = 0; i < K; ++i) { sum += C[i]; B[i] = sum; }
  } else {
    //for(S i = 0; i < K; ++i) { B[i] = sum; sum += C[i]; }
    for(S i = 0; i < K; ++i) { sum += C[i]; B[i] = sum - C[i];  }
  }
}

template <typename T, typename S>
void induceSA(const T *s, S *SA, S *C, S *B, S n, S K) {
  S *b, i, j;
  T c0, c1;
  /* compute SAl */
  if(C == B) getCounts(s, C, n, K);
  getBuckets(C, B, K, false); /* find starts of buckets */
  j = n - 1;
  b = SA + B[c1 = s[j]];
  *b++ = ((0 < j) && (s[j - 1] < c1)) ? ~j : j;
  for(i = 0; i < n; ++i) {
    j = SA[i], SA[i] = ~j;
    if(0 < j) {
      --j;
      assert(s[j] >= s[j + 1]);
      if((c0 = s[j]) != c1) { B[c1] = b - SA; b = SA + B[c1 = c0]; }
      assert(i < (b - SA));
      *b++ = ((0 < j) && (s[j - 1] < c1)) ? ~j : j;
    }
  }
  /* compute SAs */
  if(C == B) getCounts(s, C, n, K);
  getBuckets(C, B, K, true); /* find ends of buckets */
  for(i = n - 1, b = SA + B[c1 = 0]; 0 <= i; --i) {
    if(0 < (j = SA[i])) {
      --j;
      assert(s[j] <= s[j + 1]);
      if((c0 = s[j]) != c1) { B[c1] = b - SA; b = SA + B[c1 = c0]; }
      assert((b - SA) <= i);
      *--b = ((j == 0) || (s[j - 1] > c1)) ? ~j : j;
    } else {
      SA[i] = ~j;
    }
  }
}

template <typename T, typename S>
void LMSsort2(const T *s, S *SA, S *C, S *B, S *D, S n, S K) {
  S *b, i, j, t, d;
  T c0, c1;
  assert(C != B);

  /* compute SAl */
  getBuckets(C, B, K, false); /* find starts of buckets */
  j = n - 1;
  b = SA + B[c1 = s[j]];
  --j;
  t = (s[j] < c1);
  j += n;
  *b++ = (t & 1) ? ~j : j;
  for(i = 0, d = 0; i < n; ++i) {
    if(0 < (j = SA[i])) {
      if(n <= j) { d += 1; j -= n; }
      assert(s[j] >= s[j + 1]);
      if((c0 = s[j]) != c1) { B[c1] = b - SA; b = SA + B[c1 = c0]; }
      assert(i < (b - SA));
      --j;
      t = c0; t = (t << 1) | (s[j] < c1);
      if(D[t] != d) { j += n; D[t] = d; }
      *b++ = (t & 1) ? ~j : j;
      SA[i] = 0;
    } else if(j < 0) {
      SA[i] = ~j;
    }
  }
  for(i = n - 1; 0 <= i; --i) {
    if(0 < SA[i]) {
      if(SA[i] < n) {
        SA[i] += n;
        for(j = i - 1; SA[j] < n; --j) { }
        SA[j] -= n;
        i = j;
      }
    }
  }

  /* compute SAs */
  getBuckets(C, B, K, true); /* find ends of buckets */
  for(i = n - 1, d += 1, b = SA + B[c1 = 0]; 0 <= i; --i) {
    if(0 < (j = SA[i])) {
      if(n <= j) { d += 1; j -= n; }
      assert(s[j] <= s[j + 1]);
      if((c0 = s[j]) != c1) { B[c1] = b - SA; b = SA + B[c1 = c0]; }
      assert((b - SA) <= i);
      --j;
      t = c0; t = (t << 1) | (s[j] > c1);
      if(D[t] != d) { j += n; D[t] = d; }
      *--b = (t & 1) ? ~(j + 1) : j;
      SA[i] = 0;
    }
  }
}

template <typename S>
S LMSpostproc2(S *SA, S n, S m) {
  S i, j, d, name;

  /* compact all the sorted LMS substrings into the first m items of SA */
  assert(0 < n);
  for(i = 0, name = 0; (j = SA[i]) < 0; ++i) {
    j = ~j;
    if(n <= j) { name += 1; }
    SA[i] = j;
    assert((i + 1) < n);
  }
  if(i < m) {
    for(d = i, ++i;; ++i) {
      assert(i < n);
      if((j = SA[i]) < 0) {
        j = ~j;
        if(n <= j) { name += 1; }
        SA[d++] = j; SA[i] = 0;
        if(d == m) { break; }
      }
    }
  }
  if(name < m) {
    /* store the lexicographic names */
    for(i = m - 1, d = name + 1; 0 <= i; --i) {
      if(n <= (j = SA[i])) { j -= n; --d; }
      SA[m + (j >> 1)] = d;
    }
  } else {
    /* unset flags */
    for(i = 0; i < m; ++i) {
      if(n <= (j = SA[i])) { j -= n; SA[i] = j; }
    }
  }
  return name;
}

template <typename T, typename S>
S LMSpostproc1(const T *s, S *SA, S n, S m) {
  S i, j, p, q, plen, qlen, name;
  T c0, c1;
  bool diff;

  /* compact all the sorted substrings into the first m items of SA
      2*m must be not larger than n (proveable) */
  assert(0 < n);
  for(i = 0; (p = SA[i]) < 0; ++i) { SA[i] = ~p; assert((i + 1) < n); }
  if(i < m) {
    for(j = i, ++i;; ++i) {
      assert(i < n);
      if((p = SA[i]) < 0) {
        SA[j++] = ~p; SA[i] = 0;
        if(j == m) { break; }
      }
    }
  }

  /* store the length of all substrings */
  i = n - 1; j = n - 1; c0 = s[n - 1];
  do { c1 = c0; } while((0 <= --i) && ((c0 = s[i]) >= c1));
  for(; 0 <= i;) {
    do { c1 = c0; } while((0 <= --i) && ((c0 = s[i]) <= c1));
    if(0 <= i) {
      SA[m + ((i + 1) >> 1)] = j - i; j = i + 1;
      do { c1 = c0; } while((0 <= --i) && ((c0 = s[i]) >= c1));
    }
  }

  /* find the lexicographic names of all substrings */
  for(i = 0, name = 0, q = n, qlen = 0; i < m; ++i) {
    p = SA[i], plen = SA[m + (p >> 1)], diff = 1;
    if((plen == qlen) && ((q + plen) < n)) {
      for(j = 0; (j < plen) && (s[p + j] == s[q + j]); ++j) { }
      if(j == plen) { diff = 0; }
    }
    if(diff != 0) { ++name, q = p, qlen = plen; }
    SA[m + (p >> 1)] = name;
  }
  return name;
}

template<typename T, typename S>
void LMSsort1(const T *s, S *SA, S *C, S *B, S n, S K) {
  S *b, i, j;
  T c0, c1;

  /* compute SAl */
  if(C == B) { getCounts(s, C, n, K); }
  getBuckets(C, B, K, false); /* find starts of buckets */
  j = n - 1;
  b = SA + B[c1 = s[j]];
  --j;
  *b++ = (s[j] < c1) ? ~j : j;
  for(i = 0; i < n; ++i) {
    if(0 < (j = SA[i])) {
      assert(s[j] >= s[j + 1]);
      if((c0 = s[j]) != c1) { B[c1] = b - SA; b = SA + B[c1 = c0]; }
      assert(i < (b - SA));
      --j;
      *b++ = (s[j] < c1) ? ~j : j;
      SA[i] = 0;
    } else if(j < 0) {
      SA[i] = ~j;
    }
  }
  /* compute SAs */
  if(C == B) { getCounts(s, C, n, K); }
  getBuckets(C, B, K, true); /* find ends of buckets */
  for(i = n - 1, b = SA + B[c1 = 0]; 0 <= i; --i) {
    if(0 < (j = SA[i])) {
      assert(s[j] <= s[j + 1]);
      if((c0 = s[j]) != c1) { B[c1] = b - SA; b = SA + B[c1 = c0]; }
      assert((b - SA) <= i);
      --j;
      *--b = (s[j] > c1) ? ~(j + 1) : j;
      SA[i] = 0;
    }
  }
}

template <typename T, typename S>
int computeBWT(const T *s, S *SA, S *C, S *B, S n, S K)
{
  S *b, i, j;
  T c0, c1;
  int pidx = -1;
  /* compute SAl */
  if(C == B) getCounts(s, C, n, K);
  getBuckets(C, B, K, false); /* find starts of buckets */
  j = n - 1;
  b = SA + B[c1 = s[j]];
  *b++ = ((0 < j) && (s[j - 1] < c1)) ? ~j : j;
  for(i = 0; i < n; ++i) {
    if(0 < (j = SA[i])) {
      --j;
      assert(s[j] >= s[j + 1]);
      SA[i] = ~(static_cast<S>((c0 = s[j])));
      if(c0 != c1) { B[c1] = b - SA; b = SA + B[c1 = c0]; }
      assert(i < (b - SA));
      *b++ = ((0 < j) && (s[j - 1] < c1)) ? ~j : j;
    } else if(j != 0) {
      SA[i] = ~j;
    }
  }
  /* compute SAs */
  if(C == B) getCounts(s, C, n, K);
  getBuckets(C, B, K, true); /* find ends of buckets */
  for(i = n - 1, b = SA + B[c1 = 0]; 0 <= i; --i) {
    if(0 < (j = SA[i])) {
      --j;
      assert(s[j] <= s[j + 1]);
      SA[i] = (c0 = s[j]);
      if(c0 != c1) { B[c1] = b - SA; b = SA + B[c1 = c0]; }
      assert((b - SA) <= i);
      *--b = ((0 < j) && (s[j - 1] > c1)) ? ~(static_cast<S>(s[j - 1])) : j;
    } else if(j != 0) {
      SA[i] = ~j;
    } else {
      pidx = i;
    }
  }
  return pidx;
}

} //empty namespace



template <typename T, typename S>
int sais(const T *s, S *SA, S fs, S n, S K, bool isbwt) {
  assert(s); assert(SA);
  assert(fs >= 0); assert(n > 0); assert( K >= 1);
  S *C, *B, *D, *RA, name, t, *b, i, j ,m, newfs, p, q;
  T c0, c1;
  unsigned flags;
  int pidx = 0;

  if(K <= kMinBucketSize) {
    C = new S[K];
    if(K <= fs) {
      B = SA + (n + fs - K);
      flags = 1;
    } else {
      B = new S[K];
      flags = 3;
    }
  } else if (K <= fs) {
    C = SA + (n + fs - K);
    if(K <= (fs - K)) {
      B = C - K;
      flags = 0;
    } else if (K <= (kMinBucketSize * 4)) {
      B = new S[K];
      flags = 2;
    } else {
      B = C;
      flags = 8;
    }
  } else {
    C = B = new S[K];
    flags = 4 | 8;
  }

  if ((n <= kSaisLMSsort2Limit) && ((n/K) >= 2)) {
    if(flags & 1) 
      flags |= ((2*K) <= (fs - K)) ? 32 : 16;
    else if((flags == 0) && ((2*K) <= (fs - 2*K)))
      flags |= 32;
  }

  /* stage 1: sort all the LMS-substrings */
  getCounts(s, C, n, K);
  getBuckets(C, B, K, true);
  std::fill(SA, SA + n, 0);
  b = &t; i = n - 1; j = n; m = 0; c0 = s[n-1];
  do {c1 = c0; } while((--i >= 0) && ((c0 = s[i]) >= c1));
  for(; 0 <= i;) {
    do {c1 = c0; } while((--i >= 0) && ((c0 = s[i]) <= c1));
    if(i >= 0) {
      *b = j; b = SA + --B[c1]; j = i; ++m;
      do { c1 = c0; } while((--i >= 0) && ((c0 = s[i]) >= c1));
    }
  }

  if(m > 1) {
    if(flags & (16 | 32)) {
      if(flags & 16) {
        D = new S[2*K];
      } else {
        D = B - K * 2;
      }
      assert((j + 1) < n);
      ++B[s[j + 1]];
      for(i = 0, j = 0; i < K; ++i) {
        j += C[i];
        if(B[i] != j) { assert(SA[B[i]] != 0); SA[B[i]] += n; }
        D[i] = D[i + K] = 0;
      }
      LMSsort2(s, SA, C, B, D, n, K);
      name = LMSpostproc2(SA, n, m);
      if(flags & 16) delete [] D;
    } else {
      LMSsort1(s, SA, C, B, n, K);
      name = LMSpostproc1(s, SA, n, m);
    }
  } else if(m == 1) {
    *b = j + 1;
    name = 1;
  } else {
    name = 0;
  }

  /* stage 2 */
  if(name < m) {
    if(flags & 4) delete [] C;  
    if(flags & 2) delete [] B;
    newfs = (n + fs) - (m * 2);
    if((flags & (1 | 4 | 8)) == 0) {
      if((K + name) <= newfs) { newfs -= K; }
      else { flags |= 8; }
    }
    assert((n >> 1) <= (newfs + m));
    RA = SA + m + newfs;
    for(i = m + (n >> 1) - 1, j = m - 1; m <= i; --i) {
      if(SA[i] != 0) {
        RA[j--] = SA[i] - 1;
      }
    }
    if(sais(RA, SA, newfs, m, name, false, 0) != 0) {
      if(flags & 1) delete [] C;
    }

    i = n - 1; j = m - 1; c0 = s[n - 1];
    do { c1 = c0; } while((0 <= --i) && ((c0 = s[i]) >= c1));
    for(; 0 <= i;) {
      do { c1 = c0; } while((0 <= --i) && ((c0 = s[i]) <= c1));
      if(0 <= i) {
        RA[j--] = i + 1;
        do { c1 = c0; } while((0 <= --i) && ((c0 = s[i]) >= c1));
      }
    }
    for(i = 0; i < m; ++i) { SA[i] = RA[SA[i]]; }
    if(flags & 4) {
      C = B = new S[K];
    }
    if(flags & 2) {
      B = new S[K];
    }
  }

  /* stage 3 */
  if(flags & 8) getCounts(s, C, n, K);
  /* put all left-most S characters into their buckets */
  if(1 < m) {
    getBuckets(C, B, K, true); /* find ends of buckets */
    i = m - 1, j = n, p = SA[m - 1], c1 = s[p];
    do {
      q = B[c0 = c1];
      while(q < j) { SA[--j] = 0; }
      do {
        SA[--j] = p;
        if(--i < 0) { break; }
        p = SA[i];
      } while((c1 = s[p]) == c0);
    } while(0 <= i);
    while(0 < j) { SA[--j] = 0; }
  }
  if(!isbwt) { induceSA(s, SA, C, B, n, K); }
  else { pidx = computeBWT(s, SA, C, B, n, K); }
  if(flags & (1 | 4)) delete [] C;
  if(flags & 2) delete [] B;

  return pidx;
}



} //namespace sa_is

#undef isLMS

#endif
