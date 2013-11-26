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

               cout << col << endl;

               if(Operator::gsparse(i,0,col))
                  contract_op(i,col,Operator::gop(i,0,col),tmp5,grad[i-no]);

               ++col;

            }

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

      QSDindexed_contract(1.0,tmp4,shape(l,k,n,o),A.conjugate(),shape(o,i,m),0.0,E_op,shape(l,k,n,m,i));

   }

   /**
    * contract the left operator with the correct right hand operator
    */
   void contract_op(int site,int opnum,const QSDArray<2> &op,const QSDArray<5> &E_op,QSDArray<4> &G){

      enum {i,j,k,l,m,n,o};

      //first contract the left renormalized operator with the correct on-site hamiltonian term
      QSDArray<5> E;

      QSDindexed_contract(1.0,E_op,shape(l,k,n,m,i),op,shape(i,j),0.0,E,shape(k,j,n,l,m));

      //get the correct right hand side renormalized operator (complementary)
      char name[100];

      sprintf(name,"scratch/Right/site_%d/%d.mpx",site + 1,opnum);

      std::ifstream fin(name);
      boost::archive::binary_iarchive iar(fin);

      QSDArray<3> tmp3;
      iar >> tmp3;

      //and contract
      QSDindexed_contract(1.0,E,shape(k,j,n,l,m),tmp3.conjugate(),shape(l,o,m),1.0,G,shape(k,j,n,o));

   }

}
