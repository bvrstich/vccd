#include <iostream>
#include <iomanip>

using std::cout;
using std::endl;

#include "include.h"

using namespace btas;

/**
 * constructor will allocate place for L Tensor objects
 * @param L_in length of the chain
 * @param qt_in total quantumnumber of the chain
 * @param D_in the max dimension of the symmetryblocks
 */
MPS::MPS(int L_in,const Quantum &qt_in,int D_in){

   L = L_in;
   D = D_in;

   mps = new Tensor * [L];
   qt = new Quantum(qt_in);

   //function calculates the possible symmetryblocks and dimensions and allocates the memory for the tensors
   this->initialize();

}

/**
 * copy constructor will copy an MPS of length L to this
 * @param L_in length of the chain
 */
MPS::MPS(const MPS &mps_c){

   L = mps_c.gL();
   D = mps_c.gD();

   qt = new Quantum(mps_c.gqt());

   mps = new Tensor * [L];

   for(int i = 0;i < L;++i)
      mps[i] = new Tensor(mps_c[i]);

}

/**
 * destruct deallocates the memory
 */
MPS::~MPS(){

   delete qt;

   for(int i = 0;i < L;++i)
      delete mps[i];

   delete [] mps;

}

/**
 * initialize the MPS chain on total Quantum number qt
 * @param qt total quantumnumber
 * @param D maximal dimension of the quantum blocks
 */
void MPS::initialize(){

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

   qr[L-1] = Qshapes(1,*qt);
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

   mps[0] = new Tensor (qshape,dshape);

   //then the  middle ones
   for(int i = 1;i < L;++i){

      ql = qr[i - 1];
      dl = dr[i - 1];

      qshape = blitz::TinyVector<Qshapes,3>(ql,qp,-qr[i]);
      dshape = blitz::TinyVector<Dshapes,3>(dl,dp,dr[i]);

      mps[i] = new Tensor (qshape,dshape);

   }

}

/**
 * @return the length of the chain
 */
int MPS::gL() const {

   return L;

}

/**
 * @return the maximal dimension of the symmetry blocks
 */
int MPS::gD() const {

   return D;

}

/**
 * @return the Tensor on index i
 * @param i the index
 */
const Tensor &MPS::operator[](int i) const{

   return *mps[i];

}

/**
 * @return the Tensor on index i
 * @param i the index
 */
const Quantum &MPS::gqt() const{

   return *qt;

}


ostream &operator<<(ostream &output,MPS &mps_p){

   for(int i = 0;i < mps_p.gL();++i){

      output << "Tensor on block " << i << endl;
      output << endl;
      output << mps_p[i] << endl;
      output << endl;

   }

   return output;

}

/**
 * Canonicalize the MPS after random initialization using
 * @param left if true left canonicalize, if false right
 */
void MPS::canonicalize(bool left){

   for(int i = 0;i < L - 1;++i)
      mps[i]->canonicalize(left);

   mps[L - 1]->normalize();

}

/**
 * @return the Tensor on index i
 * @param i the index
 */
Tensor &MPS::operator[](int i){

   return *mps[i];

}
