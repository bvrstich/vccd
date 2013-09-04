#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstdlib>

using std::cout;
using std::endl;
using std::ofstream;
using std::ifstream;

#include "FermiQuantum.h"
namespace btas { typedef FermiQuantum Quantum; }; // Defined as default quantum number class

/**
 * simple random number generator
 */
double rgen() { return 2.0*(static_cast<double>(rand())/RAND_MAX) - 1.0; }

#include "include.h"

using namespace btas;
using namespace mps;

int main(void){

   cout.precision(10);
   srand(time(NULL));

   //lenght of the chain
   int L = 4;

   //number of particles
   int n_u = 2;
   int n_d = 2;

   int no = n_u;
   int nv = L - no;

   Ostate::construct_oplist(L);

   Qshapes<Quantum> qp;
   physical(qp);

   MPS<Quantum> A =  mps::create(L,Quantum(n_u,n_d),qp,40,rgen); 
   mps::compress(A,mps::Right,100);
   mps::compress(A,mps::Left,100);
   mps::normalize(A);

   MPS<Quantum> B = mps::create(L,Quantum(n_u,n_d),qp,40,rgen); 
   mps::compress(B,mps::Right,100);
   mps::compress(B,mps::Left,100);
   mps::normalize(B);

   DArray<2> t(L,L);
   t.generate(rgen);

   for(int i = 0;i < L;++i)
      for(int j = i + 1;j < L;++j)
         t(i,j) = t(j,i);

   t = 0.0;

   DArray<4> V(L,L,L,L);

   for(int i = 0;i < L;++i)
      for(int j = 0;j < L;++j)
         for(int k = 0;k < L;++k)
            for(int l = 0;l < L;++l){

               double value = rgen();

               V(i,j,k,l) = value;
               V(j,i,l,k) = value;
               V(k,j,i,l) = value;
               V(j,k,l,i) = value;
               V(i,l,k,j) = value;
               V(l,i,j,k) = value;
               V(k,l,i,j) = value;
               V(l,k,j,i) = value;

            }



   V = 0.0;

   for(int i = 0;i < L;++i)
      for(int j = 0;j < L;++j)
         for(int k = 0;k < L;++k){

            double value = rgen();

            V(i,i,j,k) = value;
            V(i,i,k,j) = value;
            V(j,k,i,i) = value;
            V(k,j,i,i) = value;
            V(j,i,i,k) = value;
            V(i,k,j,i) = value;
            V(i,j,k,i) = value;
            V(k,i,i,j) = value;

         }
   MPO<Quantum> qc_test = qcham_test<Quantum>(t,V);
   compress(qc_test,mps::Right,0);
   compress(qc_test,mps::Left,0);

   MPO<Quantum> qc = qcham<Quantum>(t,V);
   compress(qc,mps::Right,0);
   compress(qc,mps::Left,0);

   cout << inprod(mps::Left,A,qc,B) << endl;
   cout << inprod(mps::Left,B,qc,A) << endl;
   cout << 0.5 * inprod(mps::Left,A,qc_test,B) << endl;
   cout << 0.5 * inprod(mps::Left,B,qc_test,A) << endl;

   return 0;

}
