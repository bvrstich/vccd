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

   MPO O = nnn_hopping(L,d,0.5);

   MPS OA = gemv(O,A);

   for(int i = 0;i < L;++i){

      cout << "site " << i << endl;
      cout << endl;
      cout << OA[i].qshape() << endl;
      cout << OA[i].dshape() << endl;
      cout << endl;

   }

   clean(OA);
   compress<3>(OA,true,0,false);
   clean(OA);
   compress<3>(OA,false,0,false);

   for(int i = 0;i < L;++i){

      cout << "site " << i << endl;
      cout << endl;
      cout << OA[i].qshape() << endl;
      cout << OA[i].dshape() << endl;
      cout << endl;

   }

   return 0;

}
