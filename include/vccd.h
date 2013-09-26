#ifndef VCCD_H
#define VCCD_H

#include <iostream>
#include <iomanip>

class Ostate;
class T2_2_mpo;

using namespace btas;
using namespace mpsxx;

namespace vccd{

   template<class Q>
      void gradient(const MPO<Q> &qcham,const MPS<Q> &wccd,DArray<4> &grad);

   template<class Q>
      void gradient_new(const DArray<4> &,const MPO<Q> &qcham,double E,const MPS<Q> &wccd,const T2_2_mpo &list,DArray<4> &grad);

   template<class Q>
      double line_search(const MPO<Q> &qcham,const MPS<Q> &hf,const DArray<4> &t,const DArray<4> &dir,double guess,const std::vector<int> &cutoff);

   template<class Q>
      double line_search_func(double,const DArray<4> &t,const DArray<4> &dir,const MPO<Q> &qcham,const MPS<Q> &hf,const std::vector<int> &cutoff);

   template<class Q>
      void steepest_descent(DArray<4> &,const MPO<Q> &qcham,const MPS<Q> &hf,const std::vector<int> &cutoff);

}

#endif
