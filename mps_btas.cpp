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
   int L = 8;

   //number of particles
   int n_u = 4;
   int n_d = 4;

   int no = n_u;
   int nv = L - no;

   Ostate::construct_oplist(L);

   Qshapes<Quantum> qp;
   physical(qp);

   MPS<Quantum> A =  mps::create(L,Quantum(n_u,n_d),qp,40,rgen); 
   mps::compress(A,mps::Right,100);
   mps::compress(A,mps::Left,100);
   mps::normalize(A);

   MPS<Quantum> B =  mps::create(L,Quantum(n_u,n_d),qp,40,rgen); 
   mps::compress(B,mps::Right,100);
   mps::compress(B,mps::Left,100);
   mps::normalize(B);

   MPO<Quantum> Eop = E<Quantum>(L,2,2,1.0);

   cout << mps::dot(mps::Left,A,B) << "\t" << mps::inprod(mps::Left,A,Eop,B) << endl;
   cout << mps::dot(mps::Left,A,B) << "\t" << mps::inprod(mps::Left,B,Eop,A) << endl;

/*
   DArray<2> t(L,L);
   t.generate(rgen);

   for(int i = 0;i < L;++i)
      for(int j = i + 1;j < L;++j)
         t(i,j) = t(j,i);

   MPO<Quantum> OT = one_body<Quantum>(t);
   compress(OT,mps::Right,0);
   compress(OT,mps::Left,0);

   MPO<Quantum> OT_test = one_body_test<Quantum>(t);
   compress(OT_test,mps::Right,0);
   compress(OT_test,mps::Left,0);

   for(int i = 0;i < L;++i){

      cout << i << endl;
      cout << endl;
      cout << OT[i].qshape() << endl;
      cout << OT[i].dshape() << endl;
      cout << endl;

   }

   cout << dot(mps::Left,A,B) << "\t" << inprod(mps::Left,A,OT,B) << endl;
   cout << dot(mps::Left,A,B) << "\t" << inprod(mps::Left,A,OT_test,B) << endl;
*/
   return 0;

}
