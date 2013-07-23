#include <iostream>
#include <iomanip>

using std::cout;
using std::endl;

#include "SpinQuantum.h"
namespace btas { typedef SpinQuantum Quantum; };

#include "btas/QSDArray.h"

#include "MPS.h"

/**
 * constructor will allocate place for L btas::QSDArray<3> objects
 * @param L_in length of the chain
 */
MPS::MPS(int L_in){

   L = L_in;

   mps = new btas::QSDArray<3> * [L];

}

/**
 * copy constructor will copy an MPS of length L to this
 * @param L_in length of the chain
 */
MPS::MPS(const MPS &mps_c){

   L = mps_c.gL();
   mps = new btas::QSDArray<3> * [L];

}

/**
 * destruct deallocates the memory
 */
MPS::~MPS(){

   delete [] mps;

}

/**
 * initialize the MPS chain on total Quantum number qt
 * @param qt total quantumnumber
 * @param D maximal dimension of the quantum blocks
 */
void MPS::initialize(const btas::Quantum &qt,int D){

   //physical index
   btas::Qshapes qp;

   qp.push_back(btas::Quantum(-1));
   qp.push_back(btas::Quantum(1));

   //shape of the physical index
   btas::Dshapes dp;

   dp.push_back(1);
   dp.push_back(1);

   std::vector<btas::Qshapes> qr(L);
   std::vector<btas::Dshapes> dr(L);

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

   qr[L-1] = btas::Qshapes(1,qt);
   dr[L-1] = btas::Dshapes(1,1);

   btas::Qshapes tmpq;
   btas::Dshapes tmpd;

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

   for(int i = 0;i < L;++i){

      cout << "site " << i << endl;
      cout << endl;
      cout << qr[i] << endl;
      cout << dr[i] << endl;
      cout << endl;

   }


}

/**
 * @return the length of the chain
 */
int MPS::gL() const {

   return L;

}

/**
 * @return the Tensor on index i
 * @param i the index
 */
const btas::QSDArray<3> &MPS::operator[](int i){

   return *mps[i];

}
