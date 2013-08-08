#ifndef SPINHAMILTONIAN_H
#define SPINHAMILTONIAN_H

#include <iostream>
#include <iomanip>

template <class Q>
class Qshapes;

//some functions which initialize an MPO to a certian Hamiltonian

namespace btas{

   MPO ising(int,int,double,double);

   MPO XY(int,int,double,double);

   MPO heisenberg(int,int,double,double,double);

   MPO raise(int,int);

   MPO lower(int,int);

   MPO Sz(int,int);

   /**
    * @param d local dimension
    * @param qp Qshapes object containing the local quantumnumbers on output, input is destroyed
    */
   template<class Q>
      void physical(int d,Qshapes<Q> &qp){

         qp.clear();

         int m = -d + 1;

         while(m < d){

            qp.push_back(Q(m));

            m += 2;

         }

      }

}

#endif
