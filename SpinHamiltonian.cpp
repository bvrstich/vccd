#include <iostream>
#include <iomanip>

using std::cout;
using std::endl;
using std::ostream;

#include "include.h"

namespace btas {

   /**
    * initialize the MPO to represent a nearest-neighbour ising Hamiltonian on a lattice of size L and with coupling constant J
    * @param L length of the chain
    * @param d local dimension: i.e. defines the size of the local spins
    * @param J coupling constant
    * @param Bz magnetic field in z direction
    */
   MPO ising(int L,int d,double J,double Bz){

      double sz = 0.5 * (d - 1.0);//size of local spin

      MPO mpo(L);

      //physical indices
      Qshapes<Quantum> qp;

      physical(d,qp);

      Qshapes<Quantum> qz;
      qz.push_back(Quantum::zero());

      Qshapes<Quantum> qt;
      qt.push_back(Quantum(2));

      //incoming
      Qshapes<Quantum> qi;
      qi.push_back(Quantum::zero());//I has spin 0
      qi.push_back(Quantum(2));//S+ has spin 2

      //outgoing
      Qshapes<Quantum> qo;
      qo.push_back(Quantum::zero());//I has spin 0
      qo.push_back(Quantum(-2));//S+ has spin 2

      TVector<Qshapes<Quantum>,4> qshape = make_array(qz,qp,-qp,qo);

      //initialize the quantumnumbers of the MPO
      mpo[0].resize(Quantum::zero(),qshape);

      qshape = make_array(qi,qp,-qp,qo);

      for(int i = 1;i < L-1;++i)
         mpo[i].resize(Quantum::zero(),qshape);

      qshape = make_array(qi,qp,-qp,-qt);

      mpo[L-1].resize(Quantum::zero(),qshape);

      //first the identity
      double mz = -sz;

      for(int m = 0;m < d;++m){

         // set block elements
         DArray<4> I_op(1, 1, 1, 1);//identity
         I_op = 1.0;

         mpo[0].insert(shape(0,m,m,0),I_op);

         for(int i = 1;i < L - 1;++i){

            mpo[i].insert(shape(0,m,m,0),I_op);
            mpo[i].insert(shape(1,m,m,1),I_op);

         }

         mpo[L-1].insert(shape(1,m,m,0),I_op);

         mz += 1.0;

      }

      //then the raising operator
      mz = -sz;

      for(int m = 0; m < d - 1; ++m) {

         // set block elements
         DArray<4> Sp(1, 1, 1, 1);
         Sp = std::sqrt( (sz - mz) * (sz + mz + 1.0) );

         // insert blocks
         mpo[0].insert(shape(0,m+1,m,1),Sp);

         for(int i = 1; i < L-1; ++i) 
            mpo[i].insert(shape(0,m+1,m,1),Sp);

         mpo[L-1].insert(shape(0,m+1,m,0),Sp);

         mz  += 1.0;

      }

      return mpo;

   }

   /**
    * initialize the MPO to represent a Sz operator
    * @param L length of the chain
    * @param d local dimension: i.e. defines the size of the local spins
    */
   MPO Sz(int L,int d){

      double sz = 0.5 * (d - 1.0);//size of local spin

      MPO mpo(L);

      //physical indices
      Qshapes<Quantum> qp;
      physical(d,qp);

      Qshapes<Quantum> qz;
      qz.push_back(Quantum::zero());//Sz has spin 0

      //incoming
      Qshapes<Quantum> qi;
      qi.push_back(Quantum::zero());//Sz has spin 0
      qi.push_back(Quantum::zero());//I has spin 0

      //outgoing
      Qshapes<Quantum> qo;
      qo.push_back(Quantum::zero());//Sz has spin 0
      qo.push_back(Quantum::zero());//I has spin 0

      TVector<Qshapes<Quantum>,4> qshape = make_array(qz,qp,-qp,qo);

      //initialize the quantumnumbers of the MPO
      mpo[0].resize(Quantum::zero(),qshape);

      qshape = make_array(qi,qp,-qp,qo);

      for(int i = 1;i < L-1;++i)
         mpo[i].resize(Quantum::zero(),qshape);

      qshape = make_array(qi,qp,-qp,qz);

      mpo[L-1].resize(Quantum::zero(),qshape);

      double mz = -sz;

      for(int m = 0;m < d;++m){

         // set block elements
         DArray<4> I_op(1, 1, 1, 1);//identity
         I_op = 1.0;

         DArray<4> B_op(1, 1, 1, 1);//magnetic fieldstrength
         B_op = mz;

         mpo[0].insert(shape(0,m,m,0),I_op);
         mpo[0].insert(shape(0,m,m,1),B_op);

         for(int i = 1;i < L - 1;++i){

            mpo[i].insert(shape(0,m,m,0),I_op);
            mpo[i].insert(shape(0,m,m,1),B_op);
            mpo[i].insert(shape(1,m,m,1),I_op);

         }

         mpo[L-1].insert(shape(0,m,m,0),B_op);
         mpo[L-1].insert(shape(1,m,m,0),I_op);

         mz += 1.0;

      }

      //just a test to see if merging works
      TVector<Qshapes<Quantum>,1> qmerge;
      TVector<Dshapes,1> dmerge;

      qmerge[0] = mpo[0].qshape(3);
      dmerge[0] = mpo[0].dshape(3);

      QSTmergeInfo<1> info(qmerge,dmerge);

      QSDArray<4> tmp;
      QSTmerge(mpo[0],info,tmp);

      mpo[0] = tmp;

      for(int i = 1;i < L - 1;++i){

         //first merge the row
         qmerge[0] = mpo[i].qshape(0);
         dmerge[0] = mpo[i].dshape(0);

         info.reset(qmerge,dmerge);

         tmp.clear();

         QSTmerge(info,mpo[i],tmp);

         //then merge the column
         qmerge[0] = tmp.qshape(3);
         dmerge[0] = tmp.dshape(3);

         info.reset(qmerge,dmerge);

         mpo[i].clear();

         QSTmerge(tmp,info,mpo[i]);

      }

      //only merge row for i = L - 1
      qmerge[0] = mpo[L - 1].qshape(0);
      dmerge[0] = mpo[L - 1].dshape(0);

      info.reset(qmerge,dmerge);

      tmp.clear();

      QSTmerge(info,mpo[L - 1],tmp);

      mpo[L - 1] = tmp;

      return mpo;

   }

   /**
    * initialize the MPO to represent a spin raising operator
    * @param L length of the chain
    * @param d local dimension: i.e. defines the size of the local spins
    */
   MPO raise(int L,int d){

      double sz = 0.5 * (d - 1.0);//size of local spin

      MPO mpo(L);

      //physical indices
      Qshapes<Quantum> qp;

      physical(d,qp);

      Qshapes<Quantum> qz;
      qz.push_back(Quantum::zero());

      Qshapes<Quantum> qt;
      qt.push_back(Quantum(2));

      //incoming
      Qshapes<Quantum> qi;
      qi.push_back(Quantum::zero());//I has spin 0
      qi.push_back(Quantum(2));//S+ has spin 2

      //outgoing
      Qshapes<Quantum> qo;
      qo.push_back(Quantum::zero());//I has spin 0
      qo.push_back(Quantum(-2));//S+ has spin 2

      TVector<Qshapes<Quantum>,4> qshape = make_array(qz,qp,-qp,qo);

      //initialize the quantumnumbers of the MPO
      mpo[0].resize(Quantum::zero(),qshape);

      qshape = make_array(qi,qp,-qp,qo);

      for(int i = 1;i < L-1;++i)
         mpo[i].resize(Quantum::zero(),qshape);

      qshape = make_array(qi,qp,-qp,-qt);

      mpo[L-1].resize(Quantum::zero(),qshape);

      //first the identity
      double mz = -sz;

      for(int m = 0;m < d;++m){

         // set block elements
         DArray<4> I_op(1, 1, 1, 1);//identity
         I_op = 1.0;

         mpo[0].insert(shape(0,m,m,0),I_op);

         for(int i = 1;i < L - 1;++i){

            mpo[i].insert(shape(0,m,m,0),I_op);
            mpo[i].insert(shape(1,m,m,1),I_op);

         }

         mpo[L-1].insert(shape(1,m,m,0),I_op);

         mz += 1.0;

      }

      //then the raising operator
      mz = -sz;

      for(int m = 0; m < d - 1; ++m) {

         // set block elements
         DArray<4> Sp(1, 1, 1, 1);
         Sp = std::sqrt( (sz - mz) * (sz + mz + 1.0) );

         // insert blocks
         mpo[0].insert(shape(0,m+1,m,1),Sp);

         for(int i = 1; i < L-1; ++i) 
            mpo[i].insert(shape(0,m+1,m,1),Sp);

         mpo[L-1].insert(shape(0,m+1,m,0),Sp);

         mz  += 1.0;

      }

      return mpo;

   }


}
