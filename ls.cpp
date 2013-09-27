#include <iostream>
#include <fstream>
#include <cmath>

using std::endl;
using std::ostream;
using std::ofstream;
using std::ifstream;

#include "include.h"

using namespace btas;
using namespace mpsxx;

namespace ls {

   //constructor
   template<class Q>
      void construct(DArray<2> &Emat,DArray<2> &Nmat,const MPO<Q> &qc,const MPS<Q> &hf,const eMPS &ccd,const e_eMPS &e_emps){

         e_emps.fillE(Emat,qc,hf,ccd);
         e_emps.fillN(Nmat,ccd);

      }

   double eval(const DArray<2> &Emat,const DArray<2> &Nmat,double alpha){

      //get the faculties
      int fac[Nmat.shape(0)];

      fac[0] = 1;
      fac[1] = 1;

      for(int i = 2;i < Nmat.shape(0);++i)
         fac[i] = fac[i - 1] * i;

      double pow_alpha[Nmat.shape(0)];

      pow_alpha[0] = 1.0;
      pow_alpha[1] = alpha;

      for(int i = 2;i < Nmat.shape(0);++i)
         pow_alpha[i] = pow_alpha[i - 1] * alpha;

      double norm = 0.0;

      for(int i = 0;i < Nmat.shape(0);++i){

         norm += fac[i] * fac[i] * pow_alpha[i] * pow_alpha[i] * Nmat(i,i);

         for(int j = i + 1;j < Nmat.shape(0);++j)
            norm += 2.0 * fac[i] * fac[j] * pow_alpha[i] * pow_alpha[j] * Nmat(i,j);

      }

      double energy = 0.0;

      for(int i = 0;i < Nmat.shape(0);++i){

         energy += fac[i] * fac[i] * pow_alpha[i] * pow_alpha[i] * Emat(i,i);

         for(int j = i + 1;j < Nmat.shape(0);++j)
            energy += 2.0 * fac[i] * fac[j] * pow_alpha[i] * pow_alpha[j] * Emat(i,j);

      }

      return energy/norm;
      
   }

   template void construct<Quantum>(DArray<2> &Emat,DArray<2> &Nmat,const MPO<Quantum> &qc,const MPS<Quantum> &hf,const eMPS &emps,const e_eMPS &e_emps);

}
