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

   MPO O = n_loc(L,d,2);

   MPS OA = gemv(O,A);
   clean(OA);

   print(OA);

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
