#ifndef SPINHAMILTONIAN_H
#define SPINHAMILTONIAN_H

#include <iostream>
#include <iomanip>

template<class Q>
class Qshapes;

//some functions which initialize an MPO to a certian Hamiltonian

namespace btas{

   /**
    * @param d local dimension: number of particles on local site
    * @param qp Qshapes object containing the local quantumnumbers on output, input is destroyed
    */
   template<class Q>
      void physical(int d,Qshapes<Q> &qp){

         qp.clear();

         int n = 0;

         while(n < d){

            qp.push_back(Q(n));
            n++;

         }

      }

   MPO creator(int L,int d,int site);

   MPO annihilator(int L,int d,int site);

   MPO n_loc(int L,int d,int site);

   MPO N_tot(int L,int d,int site);

}

#endif
