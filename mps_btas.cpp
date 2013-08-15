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

#include "include.h"

using namespace btas;

int main(void){

   cout.precision(10);
   srand(time(NULL));

   //lenght of the chain
   int L = 20;

   //number of particles
   int n_u = 10;
   int n_d = 10;

   MPS A = create(L,Quantum(n_u,n_d),20);
   compress<3>(A,true,100,true);
   clean(A);

   DArray<2> t(n_u,L-n_u);
   t.generate(rgen);

   double norm = 2*Ddot(t,t);

   Dscal(1.0/sqrt(norm),t);

   MPO O = T1(t);

   compress<4>(O,true,0,false);
   clean(O);
   compress<4>(O,false,0,false);
   clean(O);

   MPS OA = gemv(O,A);

   clean(OA);
   compress<3>(OA,true,0,false);
   clean(OA);
   compress<3>(OA,false,0,false);
   clean(OA);

   return 0;

}
