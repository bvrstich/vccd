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
    * construct an MPO, i.e. a QSDArray<4> which contains a fully contracted <tccd|T H |wccd > with the T[i] missing.
    * @param no number of occupied
    * @param nv number of virtuals
    * @param E current approximation of the energy
    * @param wccd input MPS: full ccd wavefunction
    */
   MPO<Quantum> construct(int no,int nv,double E,const MPS<Quantum> &wccd){

      int L = no + nv;

      MPO<Quantum> grad(nv - no + 1);

      if(no == nv){

         //only one site to be constructed, i.e. i = no - 1
         QSDArray<5> tmp5;

         cout << endl;
         cout << no - 1 << endl;
         cout << endl;

         //get the previous operator from memory and paste the next sites onto it:
         tmp5.clear();

         //first row = 0: id coming in
         get_op(no - 2,0,wccd[no - 1],tmp5);

         int col = 0;

         while(HamOp::ostates[no - 1][col].size() == 1){

            if(Operator::gsparse(no - 1,0,col))
               contract_op(no - 1,col,Operator::gop(no - 1,0,col),tmp5,grad[0]);

            ++col;

         }

         while(HamOp::ostates[no - 1][col].size() == 2){

            if(Operator::gsparse(no - 1,0,col))
               contract_op(no - 1,col,Operator::gop(no - 1,0,col),tmp5,grad[0]);

            ++col;

         }

         int osbar = col;

         while(col < HamOp::ostates[no - 1].size()){

            if(Operator::gsparse(no - 1,0,col))
               contract_op(no - 1,col,Operator::gop(no - 1,0,col),tmp5,grad[0]);

            ++col;

         }

         for(int row = 1;row < HamOp::ostates[no - 2].size();++row){

            enum {j,k,l,m,n,o};

            cout << no - 1 << "\t" << row << "\t" << HamOp::ostates[no-2].size() << endl;

            tmp5.clear();
            get_op(no - 2,row,wccd[no - 1],tmp5);

            for(int col = 0;col < osbar;++col)
               if(Operator::gsparse(no - 1,row,col))
                  contract_op(no - 1,col,Operator::gop(no - 1,row,col),tmp5,grad[0]);

            for(int col = osbar;col < HamOp::ostates[no - 1].size();++col)
               if(Operator::gsparse(no - 1,row,col))
                  contract_op(no - 1,col,Operator::gop(no - 1,row,col),tmp5,grad[0]);

         }

         //this concludes the Hamiltonian part of the gradient, now do the -E norm
         tmp5.clear();
         get_op(no - 2,0,wccd[no - 1],tmp5);

         //read in the right identity
         QSDArray<3> tmp3;
         ro::read(mpsxx::Right,no,HamOp::ostates[no - 1].size() - 1,tmp3);

         enum {j,k,l,m,n,o};

         QSDindexed_contract(-E,tmp5,shape(k,j,n,l,m),tmp3,shape(l,o,m),1.0,grad[0],shape(k,j,n,o));

      }
      else{//no < nv

         QSDArray<5> tmp5;

         //site from no --> nv
         for(int i = no;i <= nv;++i){

            cout << endl;
            cout << i << endl;
            cout << endl;

            //get the previous operator from memory and paste the next sites onto it:
            tmp5.clear();

            //first row = 0: id coming in
            get_op(i - 1,0,wccd[i],tmp5);

            int col = 0;

            while(HamOp::ostates[i][col].size() == 1){

               if(Operator::gsparse(i,0,col))
                  contract_op(i,col,Operator::gop(i,0,col),tmp5,grad[i-no]);

               ++col;

            }

            while(HamOp::ostates[i][col].size() == 2){

               if(Operator::gsparse(i,0,col))
                  contract_op(i,col,Operator::gop(i,0,col),tmp5,grad[i-no]);

               ++col;

            }

            int osbar = col;

            while(col < HamOp::ostates[i].size()){

               if(Operator::gsparse(i,0,col))
                  contract_op(i,col,Operator::gop(i,0,col),tmp5,grad[i-no]);

               ++col;

            }

            for(int row = 1;row < HamOp::ostates[i - 1].size();++row){

               enum {j,k,l,m,n,o};

               cout << i << "\t" << row << "\t" << HamOp::ostates[i-1].size() << endl;

               tmp5.clear();
               get_op(i - 1,row,wccd[i],tmp5);

               for(int col = 0;col < osbar;++col)
                  if(Operator::gsparse(i,row,col))
                     contract_op(i,col,Operator::gop(i,row,col),tmp5,grad[i-no]);

               for(int col = osbar;col < HamOp::ostates[i].size();++col)
                  if(Operator::gsparse(i,row,col))
                     contract_op(i,col,Operator::gop(i,row,col),tmp5,grad[i-no]);

            }

            //this concludes the Hamiltonian part of the gradient, now do the -E norm
            tmp5.clear();
            get_op(i - 1,0,wccd[i],tmp5);

            //read in the right identity
            QSDArray<3> tmp3;
            ro::read(mpsxx::Right,i + 1,HamOp::ostates[i].size() - 1,tmp3);

            enum {j,k,l,m,n,o};

            QSDindexed_contract(-E,tmp5,shape(k,j,n,l,m),tmp3,shape(l,o,m),1.0,grad[i - no],shape(k,j,n,o));

         }

      }

      return grad;

   }

   /**
    * read in the operator on the previous site with number opnum, paste the next site to it:
    */
   void get_op(int site,int opnum,const QSDArray<3> &A,QSDArray<5> &E_op){

      enum {i,j,k,l,m,n,o};

      char name[100];

      sprintf(name,"scratch/Left/site_%d/%d.mpx",site,opnum);

      std::ifstream fin(name);
      boost::archive::binary_iarchive iar(fin);

      QSDArray<3> tmp3;

      iar >> tmp3;

      QSDArray<4> tmp4;

      QSDindexed_contract(1.0,tmp3,shape(o,k,m),A,shape(o,n,l),0.0,tmp4,shape(l,k,n,m));

      QSDindexed_contract(1.0,tmp4,shape(l,k,n,o),A.conjugate(),shape(o,i,m),0.0,E_op,shape(k,i,n,l,m));

   }

   /**
    * contract the left operator with the correct right hand operator
    */
   void contract_op(int site,int opnum,const QSDArray<2> &op,const QSDArray<5> &E_op,QSDArray<4> &G){

      //get the correct right hand side renormalized operator (complementary)
      char name[100];

      sprintf(name,"scratch/Right/site_%d/%d.mpx",site + 1,opnum);

      std::ifstream fin(name);
      boost::archive::binary_iarchive iar(fin);

      QSDArray<3> tmp3;
      iar >> tmp3;

      enum {i,j,k,l,m,n,o};

      QSDArray<4> tmp4;

      //and contract with the left to obtain a temporary grdient
      QSDindexed_contract(1.0,E_op,shape(k,i,n,l,m),tmp3,shape(l,o,m),0.0,tmp4,shape(k,i,n,o));

      QSDindexed_contract(1.0,tmp4,shape(k,i,n,o),op,shape(i,j),1.0,G,shape(k,j,n,o));

   }

}
