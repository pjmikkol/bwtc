/*******************************************************************
 * Burrows-Wheeler transform using the SA-IS algorithm desribed in *
 * "Two Efficient Algorithms for Linear Suffix Array Construction" *
 * by Nong, Zhang & Chan                                           *
 ******************************************************************/
#ifndef BWTC_SAIS_BWT_H_
#define BWTC_SAIS_BWT_H_

#include <cassert>

#include <algorithm>
#include <iostream>
#include <vector>

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
  for(int64 i = 0; i < n; ++i) if(isLMS(SA[i])) SA[n1++] = SA[i];

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
  for(int64 i = 0; i < n1; ++i){assert(SA1[i] != Max<S>::max); SA1[i] = s1[SA1[i]];}
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
  for(int64 i = n - 3; i >= 0; --i)
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
void InduceSAlZ(const std::vector<bool>& t, S *SA, T *s, uint32 *bkt, int64 n,
               uint32 K, bool end)
{
  GetBucketsZeroInclude(s, bkt, n, K, end);
  for(int64 i = 0; i < n; ++i) {
    S j = SA[i] - 1;
    if(j < Max<S>::max - 1 && !t[j]) {
      if(j != n - 1) SA[bkt[s[j] + 1]++] = j;
      else SA[bkt[s[j]]++] = j;
    }
  }
}

template <typename T, typename S>
void InduceSAsZ(const std::vector<bool>& t, S *SA, T *s, uint32 *bkt, int64 n,
               uint32 K, bool end)
{
  GetBucketsZeroInclude(s, bkt, n, K, end);
  for(int64 i = n - 1; i >= 0; --i) {
    S j = SA[i] - 1;
    if(j < Max<S>::max - 1 && t[j]) {
      if(j != n - 1) SA[--bkt[s[j] + 1]] = j;
      else SA[--bkt[s[j]]] = j;
    }
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
  assert(SA);
  assert(s);
  assert(s[n-1] == 0);

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
  InduceSAlZ(t, SA, s, bkt, n, K, false);
  InduceSAsZ(t, SA, s, bkt, n, K, true);

  delete [] bkt;

  /* Compact the sorted substrings into the first n1 items of SA */
  int64 n1 = 0;
  for(int64 i = 0; i < n; ++i) if(isLMS(SA[i])) SA[n1++] = SA[i];

  assert(2*n1 <= n);

  /* Find the lexicographic names of all substrings */
  std::fill(SA + n1, SA + n, Max<S>::max);
  SA[n1 + SA[0]/2] = 0;

  S name = 1, prev = Max<S>::max;
  for(int64 i = 1; i < n1; ++i) {
    S pos = SA[i];
    bool diff = false;
    for(S d = 0; d < n; ++d) {
      if(prev == Max<S>::max || s[pos+d] != s[prev+d]||t[pos+d] != t[prev+d])
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
  bkt = new uint32[K+1];
  /* Put all LMS-characters into their buckets. */
  GetBucketsZeroInclude(s, bkt, n, K, true);
  for(int64 i = 1, j = 0; i < n; ++i) 
    if(isLMS(i)) s1[j++] = i; 
  for(int64 i = 0; i < n1; ++i) SA1[i] = s1[SA1[i]];
  std::fill(SA + n1, SA + n, Max<S>::max);
  for(int64 i = n1 - 1; i >= 1; --i) {
    S j = SA[i];
    SA[i] = Max<S>::max;
    SA[--bkt[s[j] + 1]] = j;
  }
  SA[--bkt[s[SA[0]]]] = SA[0];
  InduceSAlZ(t, SA, s, bkt, n, K, false);
  InduceSAsZ(t, SA, s, bkt, n, K, true);

  delete [] bkt;
}

} //namespace sa_is

#undef isLMS

#endif
