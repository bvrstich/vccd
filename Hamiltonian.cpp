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

      MPO mpo(L);

      //physical indices
      Qshapes<Quantum> qp;

      qp.push_back(-1);
      qp.push_back(1);

      //incoming
      Qshapes<Quantum> qi;
      qi.push_back(Quantum::zero());//Sz has spin 0

      //outgoing
      Qshapes<Quantum> qo;
      qo.push_back(Quantum::zero());//Sz has spin 0

      TVector<Qshapes<Quantum>,4> qshape = make_array(qi,qp,-qp,qo);

      //initialize the quantumnumbers of the MPO
      for(int i = 0;i < L;++i)
         mpo[i].resize(Quantum::zero(),qshape);

      return mpo;

   }

}
