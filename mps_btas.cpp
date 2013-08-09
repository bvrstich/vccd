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

   //physical dimension
   int d = 2;

   //number of particles
   int N = 2;

   MPS A = create(L,d,Quantum(N),10);

   MPO O_cr = creator(L,d,2);
   MPO O_an = annihilator(L,d,1);

   MPS OA = gemv(O_an,A);
   clean(OA);

   cout << OA[1] << endl;

/*
   MPS OA2 = gemv(O_an,OA);
   clean(OA2);
   for(int i = 0;i < L;++i){

      cout << endl;
      cout << "site " << i << endl;
      cout << OA[i].qshape() << endl;
      cout << OA[i].dshape() << endl;
      cout << endl;

   }

*/
   return 0;

}
