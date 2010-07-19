#include <cassert>
#include <cstdlib>

#include <vector>

#include "../globaldefs.h"

namespace long_sequences {

namespace {
/* Utilities for finding prime numbers. Test we use here is deterministic *
 * variant of Miller-Rabin primality test. */

int mulm(int a, int b, int m)
{
  return static_cast<int64>(a)*b % m;
}

template<class T>
T powm(T b, T e, T m)
{
  if (e == 0) return 1;
  T x = powm(b, e>>1, m);
  T xx = mulm(x, x, m);
  return (e & 1) ? mulm(b, xx, m) : xx;
}

template<class T>
bool Witness(T a, T n)
{
  // TODO: Replace use of builtin for other compilers than gcc
  int t = __builtin_ffsll(n-1)-1;
  T x = powm(a, (n-1) >> t, n);
  for(int i = 0; i < t; ++i) {
    T y = mulm(x,x,n);
    if (y == 1 && x != 1 && x != n-1) return 1;
    x = y;
  }
  return x != 1;
}


} // empty namespace

bool IsPrime(int n)
{
  assert(n % 2 == 1);
  if ( n < 1373653) return (!Witness(2,n) && !Witness(3,n));
  else if ( n < 9080191 ) return (!Witness(31,n) && !Witness(73,n));
  /* n < 4,759,123,141 */
  else return (!Witness(2,n) && !Witness(7,n) && !Witness(61,n));
}

/* Finds prime which is greater or equal to given n */
int FindPrimeGeq(int n) {
  if (n%2 == 0) n += 1;
  while (!IsPrime(n)) n += 2;
  return n;
}

/* Hashtable based structure for holding the frequency information about the *
 * substrings encountered. */
class SequenceTable {
  struct seq {
    seq() : count(0), overflow(0), position(0) {}
    seq(unsigned c, int o, uint64 p) : count(c), overflow(o), position(p) {}
    
    unsigned count; /* How many times we have encountered the string */
    int overflow; /* Index of overflow sequence in over_ -vector */
    uint64 position; /* Position in actual text */
  };

 public:
  explicit SequenceTable(unsigned size, byte *data) : data_(data)
  {
    size_ = FindPrimeGeq(size);
    table_ = new seq[size_];
  }

  ~SequenceTable() {
    delete [] table_;
  }

 private:
  byte *data_; /* Source of strings stored to table */
  unsigned size_; /* Size of actual hashtable */
  seq *table_; /* hashtable */
  std::vector<seq> over_; /* overflow- lists*/
};

template <typename T>
int Border(T *source, unsigned length) {
  std::vector<int> border(length + 1);
  border[0] = -1;
  int t = -1;
  for(unsigned j = 1; j <= length; ++j) {
    while(t >= 0 && source[t] != source[j-1]) t = border[t];
    border[j] = ++t;
  }
  return border[length];
}

} //namespace long_sequences


uint64 CompressSequences(byte *from, uint64 length)
{
  using namespace long_sequences;

  assert(length > 0);
  assert(from);
  

}
