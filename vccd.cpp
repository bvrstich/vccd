#include <iostream>
#include <iomanip>

using std::cout;
using std::endl;
using std::ostream;

#include "include.h"

using namespace mps;
using namespace btas;

namespace vccd {

   /**
    * construct the gradient of the energy for ccd
    */
   template<class Q>
      void gradient(const MPO<Q> &qcham,const MPS<Q> &wccd,DArray<4> &grad){

         double norm_ccd = dot(mps::Left,wccd,wccd);
         MPS<Q> Hccd = gemv(qcham,wccd);
         compress(Hccd,mps::Right,0);
         compress(Hccd,mps::Left,0);

         double E_ccd = dot(mps::Left,wccd,Hccd)/norm_ccd;

         cout << norm_ccd << "\t" << E_ccd << endl;

         int no = grad.shape(0);//number of occupied orbitals
         int nv = grad.shape(2);//number of virtual orbitals

         int L = no + nv;

         MPO<Q> eai;
         MPO<Q> ebj;
         MPO<Q> prod;
         MPS<Q> tmp;

         for(int i = 0;i < no;++i)
            for(int j = 0;j < no;++j)
               for(int a = 0;a < nv;++a)
                  for(int b = 0;b < nv;++b){

                     grad(i,j,a,b) = 0.0;

                     eai = E<Quantum>(L,a,i,1.0);
                     ebj = E<Quantum>(L,b,j,1.0);

                     prod = gemm(eai,ebj);
                     tmp = gemv(prod,wccd);

                     grad(i,j,a,b) = dot(mps::Left,tmp,Hccd) - E_ccd * dot(mps::Left,tmp,wccd)/norm_ccd;

                  }


      }

   template void gradient<Quantum>(const MPO<Quantum> &,const MPS<Quantum> &wccd,DArray<4> &grad);

}
