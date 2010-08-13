/*******************************************************************
 * Burrows-Wheeler transform using the SA-IS algorithm desribed in *
 * "Two Efficient Algorithms for Linear Suffix Array Construction" *
 * by Nong, Zhang & Chan                                           *
 ******************************************************************/
#ifndef BWTC_SAIS_BWT_H_
#define BWTC_SAIS_BWT_H_

#include <vector>

#include "bw_transform.h"
#include "../globaldefs.h"

namespace bwtc {

class SAISBWTransform : public BWTransform {
 public:
  SAISBWTransform();
  virtual ~SAISBWTransform();

 private:
  void GetBuckets(byte *s, int *bkt, int n, int K, int cs, bool end);
  void InduceSAl(std::vector<bool> *t, int *SA, byte *s, int *bkt, int n,
                 int K,int cs, bool end);
  void InduceSAs(std::vector<bool> *t, int *SA, byte *s, int *bkt, int n,
                 int K, int cs, bool end);
  void SA_IS(byte *s, int *SA, int n, int K, int cs);

};

} //namespace bwtc

#endif
