#include <iostream>
#include <fstream>
#include <cmath>

using std::endl;
using std::ostream;
using std::ofstream;
using std::ifstream;

#include "include.h"

using namespace btas;
using namespace mps;

namespace ro{

   /**
    * construct a renormalized operator for <B|O|A>
    * @param dir left renormalized or right renormalized
    * @param A input MPS
    * @param O input MPO
    * @param B input MPS
    */
   RO construct(const MPS_DIRECTION &dir,const MPS<Quantum> &A,const MPO<Quantum> &O,const MPS<Quantum> &B){

      RO renop(O.size());

      if(dir == Left){

         enum {j,k,l,m,n,o};

         //from left to right
         QSDArray<5> loc;

         QSDindexed_contract(1.0,O[0],shape(m,n,k,o),A[0],shape(j,k,l),0.0,loc,shape(j,m,n,l,o));

         //merge 2 rows together
         TVector<Qshapes<Q>,2> qmerge;
         TVector<Dshapes,2> dmerge;

         for(int i = 0;i < 2;++i){

            qmerge[i] = loc.qshape(i);
            dmerge[i] = loc.dshape(i);

         }

         QSTmergeInfo<2> info(qmerge,dmerge);

         QSDArray<4> tmp;
         QSTmerge(info,loc,tmp);

         //this will contain the right going part
         QSDArray<3> EO;

         QSDindexed_contract(1.0,B[0].conjugate(),shape(j,k,l),tmp,shape(j,k,m,n),0.0,EO,shape(m,n,l));

         QSDArray<4> I1;
         QSDArray<4> I2;

         for(int i = 1;i < L;++i){

            I1.clear();

            QSDindexed_contract(1.0,EO,shape(j,k,l),A[i],shape(j,m,n),0.0,I1,shape(k,l,n,m));

            I2.clear();

            QSDindexed_contract(1.0,I1,shape(k,l,n,m),O[i],shape(k,j,m,o),0.0,I2,shape(l,j,n,o));

            EO.clear();

            QSDindexed_contract(1.0,I2,shape(l,j,n,o),B[i].conjugate(),shape(l,j,k),0.0,EO,shape(n,o,k));

            //bad style: if no blocks remain, return zero
            if(EO.begin() == EO.end())
               return 0.0;

         }

      }
      else{

         enum {j,k,l,m,n,o};

         //from right to left
         QSDArray<5> loc;

         QSDindexed_contract(1.0,O[L - 1],shape(j,k,l,m),A[L - 1],shape(o,l,n),0.0,loc,shape(o,j,k,n,m));

         //merge 2 columns together
         TVector<Qshapes<Q>,2> qmerge;
         TVector<Dshapes,2> dmerge;

         for(int i = 0;i < 2;++i){

            qmerge[i] = loc.qshape(3 + i);
            dmerge[i] = loc.dshape(3 + i);

         }

         QSTmergeInfo<2> info(qmerge,dmerge);

         QSDArray<4> tmp;
         QSTmerge(loc,info,tmp);

         //this will contain the left going part
         QSDArray<3> EO;
         QSDindexed_contract(1.0,tmp,shape(j,k,l,m),B[L-1].conjugate(),shape(n,l,m),0.0,EO,shape(j,k,n));

         QSDArray<4> I1;
         QSDArray<4> I2;

         for(int i = L - 2;i >= 0;--i){

            I1.clear();

            QSDindexed_contract(1.0,A[i],shape(j,k,l),EO,shape(l,m,n),0.0,I1,shape(j,k,m,n));

            I2.clear();

            QSDindexed_contract(1.0,O[i],shape(l,o,k,m),I1,shape(j,k,m,n),0.0,I2,shape(j,l,o,n));

            EO.clear();

            QSDindexed_contract(1.0,B[i].conjugate(),shape(k,o,n),I2,shape(j,l,o,n),0.0,EO,shape(j,l,k));

            //bad style: if no blocks remain, return zero
            if(EO.begin() == EO.end())
               return 0.0;

         }

      }

      return renop;

   }

}
