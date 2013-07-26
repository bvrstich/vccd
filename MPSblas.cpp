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
   MPS create(int L,const Quantum &qt,int D){ 

      //physical index
      Qshapes qp;

      qp.push_back(Quantum(-1));
      qp.push_back(Quantum(1));

      //shape of the physical index
      Dshapes dp;

      dp.push_back(1);
      dp.push_back(1);

      std::vector<Qshapes> qr(L);
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

      qr[L-1] = Qshapes(1,qt);
      dr[L-1] = Dshapes(1,1);

      Qshapes tmpq;
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
      blitz::TinyVector<Qshapes,3> qshape;
      blitz::TinyVector<Dshapes,3> dshape;

      //first 0
      Qshapes ql(1,Quantum::zero());
      Dshapes dl(ql.size(),1);

      qshape = blitz::TinyVector<Qshapes,3>(ql,qp,-qr[0]);
      dshape = blitz::TinyVector<Dshapes,3>(dl,dp,dr[0]);

      //construct an MPS
      MPS mps(L);

      mps[0].resize(Quantum::zero(),qshape,dshape,rgen);

      //then the  middle ones
      for(int i = 1;i < L;++i){

         ql = qr[i - 1];
         dl = dr[i - 1];

         qshape = blitz::TinyVector<Qshapes,3>(ql,qp,-qr[i]);
         dshape = blitz::TinyVector<Dshapes,3>(dl,dp,dr[i]);

         mps[i].resize(Quantum::zero(),qshape,dshape,rgen);

      }

      return mps;

   }

   /**
    * prints all the tensors in mps_p
    * @param mps_p input MPS
    */
   void print(const MPS &mps_p){

      for(MPS::const_iterator it = mps_p.begin();it != mps_p.end();++it){

         cout << endl;
         cout << *it << endl;
         cout << endl;

      }

   }

   /**
    * will copy mps to mps_copy
    * @param mps the MPS to be copied
    * @param mps_copy the MPS into which will be copied
    */
   void copy(const MPS &mps,MPS &mps_copy){

      mps_copy.resize(mps.size());

      for(unsigned int i = 0;i < mps.size();++i)
         QSDcopy(mps[i],mps_copy[i]);

   }

   /**
    * scale the MPS with a constant factor
    * @param alpha scalingfactor
    * @param mps the MPS to be scaled
    */
   void scal(double alpha,MPS &mps){

      for(unsigned int i = 0;i < mps.size();++i)
         QSDscal(alpha,mps[i]);

   }


   /**
    * Compress an MPS object by performing an SVD
    * @param mps is the input MPS, will be lost/overwritten by the compressed MPS
    * @param left if true left canonicalize, if false right
    * @param D if > 0   this specifies the number of states to be kept
    *          if == 0  all the states are kept
    *          if < 0 all singular values > 10^-D are kept
    */
   void compress(MPS &mps,bool left,int D){

      int L = mps.size();

      if(left) {

         DiagonalQSDArray<1> S;//singular values
         QSDArray<2> V;//V^T
         QSDArray<3> U;//U --> unitary left normalized matrix

         for(int i = 0;i < L - 1;++i){

            QSDgesvd(btas::LeftCanonical,mps[i],S,U,V,D);

            //copy unitary to mps
            QSDcopy(U,mps[i]);

            //paste S and V together
            SDdidm(S,V);

            //and multiply with mps on the next site
            U = mps[i + 1];

            //when compressing dimensions will change, so reset:
            mps[i + 1].clear();

            QSDcontract(1.0,V,shape(2),U,shape(0),0.0,mps[i + 1]);

         }

         //now normalize the last tensor

         double norm = QSDdotc(mps[L - 1],mps[L - 1]);

         QSDscal(1.0/sqrt(norm),mps[L - 1]);


      }
      else{//right

         DiagonalQSDArray<1> S;//singular values
         QSDArray<3> V;//V^T --> unitary right normalized matrix
         QSDArray<2> U;//U

         for(int i = L - 1;i > 0;--i){

            QSDgesvd(btas::LeftCanonical,mps[i],S,U,V,D);

            //copy unitary to mps
            QSDcopy(V,mps[i]);

            //paste U and S together
            SDdimd(U,S);

            //and multiply with mps on the next site
            V = mps[i - 1];

            //when compressing dimensions will change, so reset:
            mps[i - 1].clear();

            QSDcontract(1.0,V,shape(2),U,shape(0),0.0,mps[i - 1]);

         }

         //now normalize the last tensor

         double norm = QSDdotc(mps[0],mps[0]);

         QSDscal(1.0/sqrt(norm),mps[0]);

      }

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

}
