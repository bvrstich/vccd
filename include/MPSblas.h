#ifndef _BTAS_MPSBLAS_H
#define _BTAS_MPSBLAS_H 1

#include <iostream>
#include <iomanip>

using std::cout;
using std::endl;
using std::ostream;

#include "SpinQuantum.h"
namespace btas { typedef SpinQuantum Quantum; };

#include "btas/QSPARSE/QSDArray.h"

namespace btas{

   //!typedefine MPS as an std::vector< QSDArray<3> > 
   typedef std::vector< QSDArray<3> > MPS;

   //!typedefine MPO as an std::vector< QSDArray<4> > 
   typedef std::vector< QSDArray<4> > MPO;

   //some function definitions on MPS's
   MPS create(int,int,const Quantum &qt,int);

   double rgen();

   /**
    * prints all the operators in mpo_p
    * @param mpo_p input MPO
    */
   template<class MPX>
      void print(const MPX &mpx_p){

         for(int i = 0;i < mpx_p.size();++i){

            cout << "tensor on site " << i << endl;
            cout << endl;
            cout << mpx_p[i] << endl;
            cout << endl;

         }

      }

   /**
    * will copy mpx to mpx_copy
    * @param mpx the MPX to be copied
    * @param mpx_copy the MPX into which will be copied
    */
   template<class MPX>
      void copy(const MPX &mpx,MPX &mpx_copy){

         mpx_copy.resize(mpx.size());

         for(unsigned int i = 0;i < mpx.size();++i)
            QSDcopy(mpx[i],mpx_copy[i]);

      }


   /**
    * scale the MPX with a constant factor
    * @param alpha scalingfactor
    * @param mpx the MPX to be scaled
    */
   template<class MPX>
      void scal(double alpha,MPX &mpx){

         for(unsigned int i = 0;i < mpx.size();++i)
            QSDscal(alpha,mpx[i]);

      }

   /**
    * construct new MPX AB that is the sum of A + B: this is done by making a larger MPX object with larger bond dimension,
    * taking the direct sum of the individual tensors in the chain
    * @param A input MPX
    * @param B input MPX
    * @return the MPX result
    */
   template<class MPX>
      MPX add(const MPX &A,const MPX &B){

         //first check if we can sum these two:
         if(A.size() != B.size())
            BTAS_THROW(false, "Error: input MP objects do not have the same length!");

         int L = A.size();

         if(A[0].qshape(1) != B[0].qshape(1))
            BTAS_THROW(false,"Error: input MP objects do not have the same physical dimension!");

         MPX AB(L);

         QSDjoin_ledge(A[0],B[0],AB[0]);

         for(int i = 1;i < L - 1;++i)
            QSDjoin(A[i],B[i],AB[i]);

         QSDjoin_redge(A[L-1],B[L-1],AB[L-1]);

         return AB;

      }

   void compress(MPS &,bool,int);

   double dot(const MPS &,const MPS &);

   double nrm2(const MPS &);

   double dist(const MPS &,const MPS &);

   double inprod(const MPS &,const MPO &,const MPS &);

   MPS gemv(const MPO &,const MPS &);

   MPS gemv2(const MPO &,const MPS &);

   MPO gemm(const MPO &,const MPO &);

   void clean(MPS &);

}

#endif 
