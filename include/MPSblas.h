#ifndef _BTAS_MPSBLAS_H
#define _BTAS_MPSBLAS_H 1

#include <iostream>
#include <iomanip>

using std::cout;
using std::endl;
using std::ostream;

#include "FermiQuantum.h"
namespace btas { typedef FermiQuantum Quantum; };

#include "btas/QSPARSE/QSDArray.h"

namespace btas{

   //!typedefine MPS as an std::vector< QSDArray<3> > 
   typedef std::vector< QSDArray<3> > MPS;

   //!typedefine MPO as an std::vector< QSDArray<4> > 
   typedef std::vector< QSDArray<4> > MPO;

   //some function definitions on MPS's
   MPS create(int,const Quantum &qt,int);
   
   MPS HF(int,const Quantum &qt);

   double rgen();

   /**
    * prints all the operators in mpo_p
    * @param mpo_p input MPO
    */
   template<typename MPX>
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
   template<typename MPX>
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
   template<typename MPX>
      void scal(double alpha,MPX &mpx){

         int L = mpx.size();

         alpha = pow(alpha,1.0/(double)L);

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
   template<size_t N,typename MPX>
      MPX add(const MPX &A,const MPX &B){

         //first check if we can sum these two:
         if(A.size() != B.size())
            BTAS_THROW(false, "Error: input MP objects do not have the same length!");

         int L = A.size();

         MPX AB(L);

         QSDArray<N> tmp;

         IVector<N-1> left;

         for(int i = 0;i < N-1;++i)
            left[i] = i;

         //first left 
         QSDdsum(A[0],B[0],left,tmp);

         //merge the column quantumnumbers together
         TVector<Qshapes<Quantum>,1> qmerge;
         TVector<Dshapes,1> dmerge;

         qmerge[0] = tmp.qshape(N-1);
         dmerge[0] = tmp.dshape(N-1);

         QSTmergeInfo<1> info(qmerge,dmerge);

         //then merge
         QSTmerge(tmp,info,AB[0]);

         IVector<N-2> middle;

         for(int i = 1;i < N-1;++i)
            middle[i - 1] = i;

         //row and column addition in the middle of the chain
         for(int i = 1;i < L - 1;++i){

            QSDdsum(A[i],B[i],middle,AB[i]);

            //merge the row quantumnumbers together
            qmerge[0] = AB[i].qshape(0);
            dmerge[0] = AB[i].dshape(0);

            info.reset(qmerge,dmerge);

            //then merge
            QSTmerge(info,AB[i],tmp);

            //column quantumnumbers
            qmerge[0] = tmp.qshape(N-1);
            dmerge[0] = tmp.dshape(N-1);

            info.reset(qmerge,dmerge);

            //then merge
            QSTmerge(tmp,info,AB[i]);

         }

         IVector<N-1> right;

         for(int i = 0;i < N-1;++i)
            right[i] = i + 1;

         //finally the right
         tmp.clear();
         QSDdsum(A[L-1],B[L-1],right,tmp);

         //merge the row quantumnumbers together
         qmerge[0] = tmp.qshape(0);
         dmerge[0] = tmp.dshape(0);

         info.reset(qmerge,dmerge);

         //then merge
         QSTmerge(info,tmp,AB[L-1]);

         return AB;

      }

   /**
    * Compress an MP object by performing an SVD
    * @param mpx is the input MPX, will be lost/overwritten by the compressed MPX
    * @param left if true left canonicalize, if false right
    * @param D if > 0   this specifies the number of states to be kept
    *          if == 0  all the states are kept
    *          if < 0 all singular values > 10^-D are kept
    */
   template<size_t N,typename MPX>
      void compress(MPX &mpx,bool left,int D,bool norm){

         int L = mpx.size();//length of the chain

         if(left) {

            SDArray<1> S;//singular values
            QSDArray<2> V;//V^T
            QSDArray<N> U;//U --> unitary left normalized matrix

            for(int i = 0;i < L - 1;++i){

               QSDgesvd(RightArrow,mpx[i],S,U,V,D);

               //copy unitary to mpx
               QSDcopy(U,mpx[i]);

               //paste S and V together
               SDdidm(S,V);

               //and multiply with mpx on the next site
               U = mpx[i + 1];

               //when compressing dimensions will change, so reset:
               mpx[i + 1].clear();

               QSDcontract(1.0,V,shape(1),U,shape(0),0.0,mpx[i + 1]);

            }

            //now normalize the last tensor
            if(norm){

               double nrm = QSDdotc(mpx[L - 1],mpx[L - 1]);

               QSDscal(1.0/sqrt(nrm),mpx[L - 1]);

            }
            else{

               //redistribute the norm over the chain
               double nrm = sqrt(QSDdotc(mpx[L-1],mpx[L-1]));

               QSDscal(1.0/nrm,mpx[L-1]);

               scal(nrm,mpx);

            }

         }
         else{//right

            SDArray<1> S;//singular values
            QSDArray<N> V;//V^T --> unitary right normalized matrix
            QSDArray<2> U;//U

            for(int i = L - 1;i > 0;--i){

               QSDgesvd(RightArrow,mpx[i],S,U,V,D);

               //copy unitary to mpx
               QSDcopy(V,mpx[i]);

               //paste U and S together
               SDdimd(U,S);

               //and multiply with mpx on the next site
               V = mpx[i - 1];

               //when compressing dimensions will change, so reset:
               mpx[i - 1].clear();

               QSDcontract(1.0,V,shape(N-1),U,shape(0),0.0,mpx[i - 1]);

            }

            //now normalize the last tensor
            if(norm){

               double nrm = QSDdotc(mpx[0],mpx[0]);

               QSDscal(1.0/sqrt(nrm),mpx[0]);

            }
            else{

               //redistribute the norm over the chain
               double nrm = sqrt(QSDdotc(mpx[0],mpx[0]));

               QSDscal(1.0/nrm,mpx[0]);

               scal(nrm,mpx);

            }

         }

      }

   /**
    * clean up the MPX, i.e. make sure the right quantumblocks are connected, remove unnecessary quantumnumbers and blocks
    * @param mpx input MPX, will be changed 'cleaned' on exit
    */
   template<typename MPX>
      void clean(MPX &mpx){

         int nlegs = mpx[0].qshape().size();

         Dshapes dr;

         int i = 0;

         //from left to right
         for(int i = 0;i < mpx.size() - 1;++i){

            dr = mpx[i].dshape()[nlegs - 1];

            std::vector<Quantum> qrem;

            for(int j = 0;j < dr.size();++j)
               if(dr[j] == 0)
                  qrem.push_back(mpx[i].qshape()[nlegs - 1][j]);//what is the quantumnumber with 0 dimension?

            if(qrem.size() != 0){

               //remove the zero blocks from site i
               for(int j = 0;j < qrem.size();++j){

                  //find the index corresponding to quantumnumber qrem[j]
                  Qshapes<Quantum> qr = mpx[i].qshape()[nlegs - 1];

                  for(int k = 0;k < qr.size();++k)
                     if(qr[k] == qrem[j])
                        mpx[i].erase(nlegs - 1,k);

               }

               for(int j = 0;j < qrem.size();++j){

                  //remove the corresponding blocks on the 0 leg of the next site
                  Qshapes<Quantum> ql = mpx[i + 1].qshape()[0];

                  for(int k = 0;k < ql.size();++k)
                     if(ql[k] == -qrem[j])
                        mpx[i + 1].erase(0,k);

               }

            }

         }

         //and back from right to left
         for(int i = mpx.size() - 1;i > 0;--i){

            dr = mpx[i].dshape()[0];//actually dl now

            std::vector<Quantum> qrem;

            for(int j = 0;j < dr.size();++j)
               if(dr[j] == 0)
                  qrem.push_back(mpx[i].qshape()[0][j]);//what is the quantumnumber with 0 dimension?

            if(qrem.size() != 0){

               //remove the zero blocks from site i
               for(int j = 0;j < qrem.size();++j){

                  //find the index corresponding to quantumnumber qrem[j]
                  Qshapes<Quantum> qr = mpx[i].qshape()[0];

                  for(int k = 0;k < qr.size();++k)
                     if(qr[k] == qrem[j])
                        mpx[i].erase(0,k);

               }

               for(int j = 0;j < qrem.size();++j){

                  //remove the corresponding blocks on the (nlegs-1) leg of the previous site
                  Qshapes<Quantum> ql = mpx[i - 1].qshape()[nlegs-1];

                  for(int k = 0;k < ql.size();++k)
                     if(ql[k] == -qrem[j])
                        mpx[i - 1].erase(nlegs-1,k);

               }

            }

         }

      }

   double dot(const MPS &,const MPS &);

   double nrm2(const MPS &);

   double dist(const MPS &,const MPS &);

   double inprod(const MPS &,const MPO &,const MPS &);

   MPS gemv(const MPO &,const MPS &);

   MPO gemm(const MPO &,const MPO &);

}

#endif 
