#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstdlib>

using std::cout;
using std::endl;
using std::ofstream;
using std::ifstream;

#include "SpinQuantum.h"
namespace btas { typedef SpinQuantum Quantum; }; // Defined as default quantum number class

#include "include.h"

using namespace btas;

int main(void){

   cout.precision(10);
   srand(time(NULL));

   //lenght of the chain
   int L = 4;

   //physical dimension
   int d = 2;

   MPS A = create(L,d,Quantum::zero(),20);

   MPO O = raise(L,d);

   MPS OA = gemv(O,A);

   for(int i = 0;i < L;++i){

      cout << endl;
      cout << OA[i].qshape() << endl;
      cout << OA[i].dshape() << endl;
      cout << endl;

   }

   clean(OA);

   cout << "CLEANED" << endl;

   for(int i = 0;i < L;++i){

      cout << endl;
      cout << OA[i].qshape() << endl;
      cout << OA[i].dshape() << endl;
      cout << endl;

   }


   return 0;

}
