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

      double Sz = 0.5 * (d - 1.0);//size of local spin

      MPO mpo(L);

      //physical indices
      Qshapes<Quantum> qp;

      qp.push_back(-1);
      qp.push_back(1);

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

      double mz = -Sz;

      for(int m = 0;m < d;++m){

         // set block elements
         DArray<4> I_op(1, 1, 1, 1);//identity
         I_op = 1.0;

         DArray<4> B_op(1, 1, 1, 1);//magnetic fieldstrength
         B_op = Bz * mz;

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

      return mpo;

   }

}
