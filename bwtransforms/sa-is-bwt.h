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
  void GetBuckets(byte *s, uint32 *bkt, uint32 n, uint32 K, int cs, bool end);
  void InduceSAl(const std::vector<bool>& t, uint32 *SA, byte *s, uint32 *bkt,
                 uint32 n, uint32 K, int cs, bool end);
  void InduceSAs(const std::vector<bool>& t, uint32 *SA, byte *s, uint32 *bkt,
                 uint32 n, uint32 K, int cs, bool end);
  void SA_IS(byte *s, uint32 *SA, uint32 n, uint32 K, int cs);

};

} //namespace bwtc

#endif
