#ifndef SPINHAMILTONIAN_H
#define SPINHAMILTONIAN_H

#include <iostream>
#include <iomanip>

template<class Q>
class Qshapes;

//some functions which initialize an MPO to a certian Hamiltonian

namespace btas{

   /**
    * @param qp Qshapes object containing the local quantumnumbers on output, input is destroyed
    */
   template<class Q>
      void physical(Qshapes<Q> &qp){

         qp.clear();

         qp.push_back(Quantum(0,0));
         qp.push_back(Quantum(0,1));
         qp.push_back(Quantum(1,0));
         qp.push_back(Quantum(1,1));

      }

   MPO creator(int L,int site,int spin);

   MPO annihilator(int L,int site,int spin);

   MPO n_loc(int L,int site);
/*
   MPO N_tot(int L,int d);
   
   MPO hopping(int L,int d);

   MPO nnn_hopping(int L,int d,double);

   MPO one_body(int L,int d,const DArray<2> &);
*/
}

#endif
