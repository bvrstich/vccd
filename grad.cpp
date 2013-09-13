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

namespace grad {

   /**
    * construct an MPO, i.e. a QSDArray<4> which contains a fully contracted <B|O|A> with the O[i] missing.
    * @param rol MPS of length L - 1 containing the left renormalized blocks
    * @param ror MPS of length L - 1 containing the right renormalized blocks
    * @param A input MPS
    * @param B input MPS
    */
   MPO<Quantum> construct(const MPS<Quantum> &rol,const MPS<Quantum> &ror,const MPS<Quantum> &A,const MPS<Quantum> &B){

      int L = A.size();

      MPO<Quantum> grad(L);

      enum {j,k,l,m,n,o};

      //first the left one
      QSDArray<4> I1;
      QSDArray<5> I2;

      QSDindexed_contract(1.0,A[0],shape(j,k,l),ror[0],shape(l,m,n),0.0,I1,shape(j,k,m,n));

      QSDindexed_contract(1.0,B[0].conjugate(),shape(o,l,n),I1,shape(j,k,m,n),0.0,I2,shape(j,o,l,k,m));

      //merge j and o together
      TVector<Qshapes<Quantum>,2> qmerge;
      TVector<Dshapes,2> dmerge;

      for(int i = 0;i < 2;++i){

         qmerge[i] = I2.qshape(i);
         dmerge[i] = I2.dshape(i);

      }

      QSTmergeInfo<2> info(qmerge,dmerge);
      QSTmerge(info,I2,grad[0]);

      //middle ones
      for(int i = 1;i < L-1;++i){

         I1.clear();
         QSDindexed_contract(1.0,rol[i-1],shape(j,k,l),A[i],shape(j,m,n),0.0,I1,shape(n,m,k,l));

         I2.clear();
         QSDindexed_contract(1.0,I1,shape(n,m,k,l),B[i].conjugate(),shape(l,j,o),0.0,I2,shape(n,m,k,j,o));

         QSDindexed_contract(1.0,I2,shape(n,m,k,j,o),ror[i],shape(n,l,o),0.0,grad[i],shape(k,j,m,l));

      }

      //last one
      I1.clear();
      QSDindexed_contract(1.0,rol[L-2],shape(j,k,l),A[L-1],shape(j,m,n),0.0,I1,shape(n,m,k,l));

      I2.clear();
      QSDindexed_contract(1.0,I1,shape(n,m,k,l),B[L-1].conjugate(),shape(l,j,o),0.0,I2,shape(k,j,m,o,n));

      //merge 2 columns together
      for(int i = 0;i < 2;++i){

         qmerge[i] = I2.qshape(3 + i);
         dmerge[i] = I2.dshape(3 + i);

      }

      info.reset(qmerge,dmerge);
      QSTmerge(I2,info,grad[L - 1]);

      return grad;

   }

   /**
    * check if the renormalized operators are correctly constructed
    */
   void check(const MPO<Quantum> &grad,const MPO<Quantum> &O){

      for(int i = 0;i < O.size();++i)
         cout << QSDdotc(grad[i].conjugate(),O[i]) << endl;

   }

}
