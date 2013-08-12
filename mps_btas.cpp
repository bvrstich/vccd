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

   //physical dimension
   int d = 2;

   //number of particles
   int N = 11;

   MPS A = create(L,d,Quantum(N),20);
   compress<3>(A,true,100,true);
   clean(A);

   MPS B = create(L,d,Quantum(N),20);
   compress<3>(B,true,100,true);
   clean(B);

   DArray<2> T(shape(L,L));
   T.generate(rgen);


   MPO O = one_body(L,d,T);

   clean(O);
   compress<4>(O,true,0,false);
   clean(O);
   compress<4>(O,false,0,false);
   clean(O);

   MPS OA = gemv(O,A);

   return 0;

}
