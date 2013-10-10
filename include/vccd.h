#ifndef VCCD_H
#define VCCD_H

#include <iostream>
#include <iomanip>

class T2_2_mpo;

using namespace btas;
using namespace mpsxx;

namespace vccd{

   template<class Q>
      void gradient(const DArray<4> &,const MPO<Q> &qcham,double E,const MPS<Q> &,const MPS<Q> &,const T2_2_mpo &list,DArray<4> &grad,bool merged);

   template<class Q>
      void solve(DArray<4> &,const MPO<Q> &qcham,const MPS<Q> &,const std::vector<double> &,int,double);

   template<class Q>
      double line_search(const MPO<Q> &qc,const MPS<Q> &hf,const DArray<4> &t,const DArray<4> &dir,double guess,int D);

}

#endif
