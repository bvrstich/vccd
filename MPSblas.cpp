#include <iostream>
#include <iomanip>

using std::cout;
using std::endl;
using std::ostream;

#include "include.h"

namespace btas {

   /**
    * create an MPS chain of length L initialized randomly on total Quantum number qt
    * @param L length of the chain
    * @param qt total quantumnumber
    * @param D maximal dimension of the quantum blocks
    * @return the MPS chain randomly filled and with correct quantumnumbers and dimensions
    */
   MPS create(int L,int d,const Quantum &qt,int D){ 

      //physical index
      Qshapes<Quantum> qp;

      physical(d,qp);

      //shape of the physical index
      Dshapes dp(qp.size(),1);

      std::vector< Qshapes<Quantum> > qr(L);
      std::vector<Dshapes> dr(L);

      qr[0] = qp;
      dr[0] = dp;

      for(int i = 1;i < L - 1;++i){

         qr[i] = qr[i-1] * qp;

         for(unsigned int j = 0;j < dr[i - 1].size();++j)
            for(unsigned int k = 0;k < dp.size();++k)
               dr[i].push_back(dr[i-1][j]*dp[k]);

         int j = 0;

         while(j < qr[i].size()){

            int k = j + 1;

            while(k < qr[i].size()){

               //this removes redundant
               if( qr[i][k] == qr[i][j] ){

                  //erase the redundant quantumnumber
                  qr[i].erase(qr[i].begin() + k);

                  //add the dimension to the right block
                  dr[i][j] += dr[i][k];

                  //erase the redundant dimension
                  dr[i].erase(dr[i].begin() + k);

               }
               else
                  ++k;

               //if dimension is too large, set to D
               if(dr[i][j] > D)
                  dr[i][j] = D;

            }

            ++j;

         }

      }

      qr[L-1] = Qshapes<Quantum>(1,qt);
      dr[L-1] = Dshapes(1,1);

      Qshapes<Quantum> tmpq;
      Dshapes tmpd;

      int i = L-2;

      for(int i = L - 2;i > 0;--i){

         tmpq = qr[i+1] * qp;
         tmpd.clear();

         for(unsigned int j = 0;j < dr[i+1].size();++j)
            for(unsigned int k = 0;k < dp.size();++k)
               tmpd.push_back(dr[i+1][j]*dp[k]);

         int j = 0;

         while(j < tmpq.size()){

            int k = j + 1;

            while(k < tmpq.size()){

               //this removes redundant
               if( tmpq[k] == tmpq[j] ){

                  //erase the redundant quantumnumber
                  tmpq.erase(tmpq.begin() + k);

                  //add the dimension to the right block
                  tmpd[j] += tmpd[k];

                  //erase the redundant dimension
                  tmpd.erase(tmpd.begin() + k);

               }
               else
                  ++k;

               //if dimension is too large, set to D
               if(tmpd[j] > D)
                  tmpd[j] = D;

            }

            ++j;

         }

         //remove irrelevant quantum blocks from below
         j = 0;

         int flag_below = 0;

         while(qr[i][j] < tmpq[0]){

            qr[i].erase(qr[i].begin() + j);
            dr[i].erase(dr[i].begin() + j);

            flag_below = 1;

         }

         //remove irrelevant quantum blocks from above
         j = qr[i].size() - 1;

         int flag_above = 0;

         while(qr[i][j] > tmpq[tmpq.size() - 1]){

            qr[i].erase(qr[i].begin() + j);
            dr[i].erase(dr[i].begin() + j);

            flag_above = 1;

            --j;

         }

         //now replace the dimensions
         for(unsigned int k = 0;k < qr[i].size();++k){

            //is there a quantumnumber in tmpq equal to qr[i][k]?
            for(unsigned int l = 0;l < tmpq.size();++l){

               //if there is, take the smallest dimension
               if(qr[i][k] == tmpq[l]){

                  if(dr[i][k] > tmpd[l])
                     dr[i][k] = tmpd[l];

               }

            }

         }

      }

      //now allocate the tensors!
      TVector<Qshapes<Quantum>,3> qshape;
      TVector<Dshapes,3> dshape;

      //first 0
      Qshapes<Quantum> ql(1,Quantum::zero());
      Dshapes dl(ql.size(),1);

      qshape = make_array(ql,qp,-qr[0]);
      dshape = make_array(dl,dp,dr[0]);

      //construct an MPS
      MPS mps(L);

      mps[0].resize(Quantum::zero(),qshape,dshape);
      mps[0].generate(rgen);

      //then the  middle ones
      for(int i = 1;i < L;++i){

         ql = qr[i - 1];
         dl = dr[i - 1];

         qshape = make_array(ql,qp,-qr[i]);
         dshape = make_array(dl,dp,dr[i]);

         mps[i].resize(Quantum::zero(),qshape,dshape);
         mps[i].generate(rgen);

      }

      return mps;

   }

   /**
    * simple random number generator
    */
   double rgen() { 

      return 2.0*(static_cast<double>(rand())/RAND_MAX) - 1.0; 

   }

   /**
    * the contraction of two MPS's
    * @return the overlap of two MPS objects
    * @param mps_X input MPS
    * @param mps_Y input MPS
    */
   double dot(const MPS &mps_X,const MPS &mps_Y){

      if(mps_X.size() != mps_Y.size())
         cout << "Error: input MPS objects do not have the same length!" << endl;

      if(mps_X[mps_X.size()-1].qshape(2) != mps_Y[mps_Y.size()-1].qshape(2))
         cout << "Error: input MPS objects do not have the same total quantumnumbers!" << endl;

      //going from left to right, this will store the already contracted part
      QSDArray<2> E;

      QSDcontract(1.0,mps_X[0],shape(0,1),mps_Y[0].conjugate(),shape(0,1),0.0,E);

      //this will contain an intermediate
      QSDArray<3> I;

      for(unsigned int i = 1;i < mps_X.size();++i){

         //construct intermediate, i.e. past mps_X to E
         QSDcontract(1.0,E,shape(0),mps_X[i],shape(0),0.0,I);

         //clear structure of E
         E.clear();

         //construct E for site i by contracting I with mps_Y
         QSDcontract(1.0,I,shape(0,1),mps_Y[i].conjugate(),shape(0,1),0.0,E);

         I.clear();

      }

      return (*(E.begin()->second))(0,0);

   }

   /**
    * @return the norm of the state
    */
   double nrm2(const MPS &mps){

      return dot(mps,mps);

   }

   /**
    * @return the distance between 2 mps's ||X - Y||_2
    */
   double dist(const MPS  &mps_X,const MPS &mps_Y){

      return nrm2(mps_X) + nrm2(mps_Y) - 2.0 * dot(mps_X,mps_Y);

   }

   /**
    * @param A input MPS
    * @param O input MPO
    * @param B input MPS
    * @return the number containing < A | O | B >
    */
   double inprod(const MPS &A,const MPO &O,const MPS &B){

      //first check if we can sum these two:
      if(A.size() != B.size() || A.size() != O.size())
         BTAS_THROW(false, "Error: input objects do not have the same length!");

      int L = A.size();

      //from left to right
      QSDArray<5> loc;

      QSDcontract(1.0,O[0],shape(2),A[0],shape(1),0.0,loc);

      //merge 2 columns together
      TVector<Qshapes<Quantum>,2> qmerge;
      TVector<Dshapes,2> dmerge;

      for(int i = 0;i < 2;++i){

         qmerge[i] = loc.qshape(3 + i);
         dmerge[i] = loc.dshape(3 + i);

      }

      QSTmergeInfo<2> info(qmerge,dmerge);

      QSDArray<4> tmp;
      QSTmerge(loc,info,tmp);

      //this will contain the right going part
      QSDArray<3> EO;

      QSDcontract(1.0,B[0].conjugate(),shape(0,1),tmp,shape(0,1),0.0,EO);

      QSDArray<4> I1;
      QSDArray<4> I2;

      for(int i = 1;i < L;++i){

         enum {j,k,l,m,n,o};

         I1.clear();

         QSDindexed_contract(1.0,EO,shape(j,k,l),A[i],shape(l,m,n),0.0,I1,shape(j,k,m,n));

         I2.clear();

         QSDindexed_contract(1.0,I1,shape(j,k,m,n),O[i],shape(k,o,m,l),0.0,I2,shape(j,o,l,n));

         EO.clear();

         QSDindexed_contract(1.0,I2,shape(j,o,l,n),B[i].conjugate(),shape(j,o,k),0.0,EO,shape(k,l,n));

      }

      return (*(EO.begin()->second))(0,0,0);

   }

   /**
    * MPO/S equivalent of a matrix vector multiplication. Let an MPO act on an MPS and return the new MPS
    * @param O input MPO
    * @param A input MPS
    * @return the new MPS object created by the multiplication
    */
   MPS gemv(const MPO &O,const MPS &A){

      //first check if we can sum these two:
      if(O.size() != A.size())
         BTAS_THROW(false, "Error: input objects do not have the same length!");

      int L = A.size();

      MPS mps(L);

      enum {j,k,l,m,n,o};

      QSDArray<5> tmp;
      QSDArray<4> mrows;

      for(int i = 0;i < L;++i){

         //clear the tmp object first
         tmp.clear();

         QSDindexed_contract(1.0,O[i],shape(j,k,l,m),A[i],shape(n,l,o),0.0,tmp,shape(n,j,k,m,o));

         //merge 2 rows together
         TVector<Qshapes<Quantum>,2> qmerge;
         TVector<Dshapes,2> dmerge;

         for(int r = 0;r < 2;++r){

            qmerge[r] = tmp.qshape(r);
            dmerge[r] = tmp.dshape(r);

         }

         QSTmergeInfo<2> info(qmerge,dmerge);

         //clear the mrows object first
         mrows.clear();

         //then merge
         QSTmerge(info,tmp,mrows);

         //merge 2 columns together
         for(int r = 2;r < 4;++r){

            qmerge[r - 2] = mrows.qshape(r);
            dmerge[r - 2] = mrows.dshape(r);

         }

         info.reset(qmerge,dmerge);

         QSTmerge(mrows,info,mps[i]);

      }

      return mps;

   }

   /**
    * MPO equivalent of a matrix matrix multiplication. MPO action on MPO gives new MPO: O1-O2|MPS>
    * @param O1 input MPO
    * @param O2 input MPO
    * @return the new MPO object created by the multiplication
    */
   MPO gemm(const MPO &O1,const MPO &O2){

      //first check if we can sum these two:
      if(O1.size() != O2.size())
         BTAS_THROW(false, "Error: input objects do not have the same length!");

      int L = O1.size();

      MPO mpo(L);

      enum {j,k,l,m,n,o,p};

      QSDArray<6> tmp;
      QSDArray<5> mrows;

      for(int i = 0;i < L;++i){

         //clear the tmp object first
         tmp.clear();

         QSDindexed_contract(1.0,O1[i],shape(n,o,k,p),O2[i],shape(j,k,l,m),0.0,tmp,shape(n,j,o,l,p,m));

         //merge 2 rows together
         TVector<Qshapes<Quantum>,2> qmerge;
         TVector<Dshapes,2> dmerge;

         for(int r = 0;r < 2;++r){

            qmerge[r] = tmp.qshape(r);
            dmerge[r] = tmp.dshape(r);

         }

         QSTmergeInfo<2> info(qmerge,dmerge);

         //clear the mrows object first
         mrows.clear();

         //then merge
         QSTmerge(info,tmp,mrows);

         //merge 2 columns together
         for(int r = 3;r < 5;++r){

            qmerge[r - 3] = mrows.qshape(r);
            dmerge[r - 3] = mrows.dshape(r);

         }

         info.reset(qmerge,dmerge);
         QSTmerge(mrows,info,mpo[i]);

      }

      return mpo;

   }

   MPS add(const MPS &A,const MPS& B){

      //first check if we can sum these two:
      if(A.size() != B.size())
         BTAS_THROW(false, "Error: input MP objects do not have the same length!");

      int L = A.size();

      if(A[0].qshape(1) != B[0].qshape(1))
         BTAS_THROW(false,"Error: input MP objects do not have the same physical dimension!");

      MPS AB(L);

      QSDArray<3> tmp;

      //first left 
      QSDdsum(A[0],B[0],shape(0,1),tmp);

      //merge the column quantumnumbers together
      TVector<Qshapes<Quantum>,1> qmerge;
      TVector<Dshapes,1> dmerge;

      qmerge[0] = tmp.qshape(2);
      dmerge[0] = tmp.dshape(2);

      QSTmergeInfo<1> info(qmerge,dmerge);

      //then merge
      QSTmerge(tmp,info,AB[0]);

      //row and column addition in the middle of the chain
      for(int i = 1;i < L - 1;++i){

         QSDdsum(A[i],B[i],shape(1),AB[i]);

         //merge the row quantumnumbers together
         qmerge[0] = AB[i].qshape(0);
         dmerge[0] = AB[i].dshape(0);

         info.reset(qmerge,dmerge);

         //then merge
         QSTmerge(info,AB[i],tmp);

         //column quantumnumbers
         qmerge[0] = tmp.qshape(2);
         dmerge[0] = tmp.dshape(2);

         info.reset(qmerge,dmerge);

         //then merge
         QSTmerge(tmp,info,AB[i]);

      }

      //finally the right
      tmp.clear();
      QSDdsum(A[L-1],B[L-1],shape(1,2),tmp);

      //merge the row quantumnumbers together
      qmerge[0] = tmp.qshape(0);
      dmerge[0] = tmp.dshape(0);

      info.reset(qmerge,dmerge);

      //then merge
      QSTmerge(info,tmp,AB[L-1]);

      return AB;

   }

}
