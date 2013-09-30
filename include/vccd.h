#ifndef VCCD_H
#define VCCD_H

#include <iostream>
#include <iomanip>

class Ostate;
class eMPS;
class T2_2_mpo;

using namespace btas;
using namespace mpsxx;

namespace vccd{

   template<class Q>
      void gradient(const DArray<4> &,const MPO<Q> &qcham,double E,const MPS<Q> &wccd,const T2_2_mpo &list,DArray<4> &grad);

   template<class Q>
      double line_search(double,double,const MPO<Q> &qcham,const MPS<Q> &,const eMPS &,const MPO<Q> &dir,double guess);

   template<class Q>
      void steepest_descent(DArray<4> &,const MPO<Q> &qcham,const MPS<Q> &hf,const std::vector<int> &cutoff);

}

#endif
