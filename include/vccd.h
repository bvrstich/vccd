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
      void gradient(const DArray<4> &,const MPO<Q> &qcham,double E,const MPS<Q> &,const MPS<Q> &,const T2_2_mpo &list,DArray<4> &grad,bool merged);

   template<class Q>
      void old_gradient(const MPO<Q> &qcham,const MPS<Q> &wccd,DArray<4> &grad,const std::vector<int> &);

   template<class Q>
      void num_gradient(double,const MPO<Q> &qcham,const MPS<Q> &hf,const DArray<4> &t,DArray<4> &grad,const std::vector<int> &);

   template<class Q>
      double line_search(double,double,const MPO<Q> &qcham,const MPS<Q> &,const eMPS &,const MPO<Q> &dir,double guess);

   template<class Q>
      void steepest_descent(DArray<4> &,const MPO<Q> &qcham,const MPS<Q> &hf,const std::vector<int> &cutoff);

   template<class Q>
      void conjugate_gradient(DArray<4> &,const MPO<Q> &qcham,const MPS<Q> &hf,const std::vector<int> &cutoff);

   template<class Q>
      void conjugate_gradient_tmp(DArray<4> &,const MPO<Q> &qcham,const MPS<Q> &hf,const std::vector<int> &cutoff);

   template<class Q>
      double old_line_search(const MPO<Q> &qc,const MPS<Q> &hf,const DArray<4> &t,const DArray<4> &dir,double guess,const std::vector<int> &cutoff);

   template<class Q>
      double line_search_func(double a,const DArray<4> &t,const DArray<4> &dir,const MPO<Q> &qc,const MPS<Q> &hf,const std::vector<int> &cutoff);

}

#endif
