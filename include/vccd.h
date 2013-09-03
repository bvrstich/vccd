#ifndef VCCD_H
#define VCCD_H

#include <iostream>
#include <iomanip>

class Ostate;

using namespace btas;
using namespace mps;

namespace vccd{

   template<class Q>
      void gradient(const MPO<Q> &qcham,const MPS<Q> &wccd,DArray<4> &grad);

}

#endif
