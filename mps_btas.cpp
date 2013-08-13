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
   int L = 4;

   //number of particles
   int n_u = 2;
   int n_d = 2;

   MPS A = create(L,Quantum(n_u,n_d),20);
   compress<3>(A,true,100,true);
   clean(A);

   MPS B = create(L,Quantum(n_u,n_d),20);
   compress<3>(B,true,100,true);
   clean(B);

   MPO O = hubbard(L,1.0);

   MPS OA = gemv(O,A);

   cout << dot(B,OA) << endl;

   clean(O);
   compress<4>(O,true,0,false);
   clean(O);
   compress<4>(O,false,0,false);
   clean(O);

   OA = gemv(O,A);

   cout << dot(B,OA) << endl;
   print(O);

   return 0;

}
