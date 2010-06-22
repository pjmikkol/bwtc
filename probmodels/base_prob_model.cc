#include <algorithm> // for std::fill
#include <iostream>

#include "../globaldefs.h"
#include "base_prob_model.h"

namespace bwtc {

ProbabilityModel* GiveProbabilityModel(char choice) {
  switch(choice) {
    case 'm':
      if( verbosity > 1)
        std::clog << "Remembering 8 previous bits\n";
      return new SimpleMarkov<byte>();
    case 'M':
      if( verbosity > 1)
        std::clog << "Remembering 16 previous bits\n";
      return new SimpleMarkov<unsigned short int>();      
    case 'n':
    default:
      if( verbosity > 1)
        std::clog << "Remembering single previous bit\n";
      return new ProbabilityModel();
  }
}


/*************************************************************************
 * SimpleMarkov: An example of how to integrate new probability model    *
 * to program.                                                           *
 *                                                                       * 
 * Simple template-based probability-model which remembers               *
 * 8*sizeof(Integer) previous bits. Consumes huge amount of memory       *
 * 2^(8*sizeof(Integer)) so practically this is  usable  only with bytes *
 * and short integers.                                                   *
 *************************************************************************/
template <typename UnsignedInt>
SimpleMarkov<UnsignedInt>::SimpleMarkov() :
    prev_(static_cast<UnsignedInt>(0)), history_(NULL)
{
  uint64 size = (static_cast<uint64>(1) << 8*sizeof(UnsignedInt)) - 1;
  history_ = new bool[size];
  std::fill(history_, history_ + size, true);
}

template <typename UnsignedInt>
SimpleMarkov<UnsignedInt>::~SimpleMarkov() {
  delete [] history_;
}

template <typename UnsignedInt>
void SimpleMarkov<UnsignedInt>::Update(bool bit) {
  history_[prev_] = bit;
  prev_ <<= 1;
  prev_ |= (bit)? 1 : 0;
}

template <typename UnsignedInt>
Probability SimpleMarkov<UnsignedInt>::ProbabilityOfOne() {
  if ( history_[prev_] ) return kProbabilityScale - 1;
  else return 1;
}

template <typename UnsignedInt>
void SimpleMarkov<UnsignedInt>::ResetModel() {
  /* Seems to work better when not resetting the model for different
   * contexts. */
  /* Resetting the model would be like.. :
     uint64 size = (static_cast<uint64>(1) << 8*sizeof(UnsignedInt)) - 1;
      std::fill(history_, history_ + size, true);*/
}

} //namespace bwtc
