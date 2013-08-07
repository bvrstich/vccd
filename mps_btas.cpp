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
   int L = 10;

   //physical dimension
   int d = 2;

   MPS A = create(L,d,Quantum::zero(),20);

   MPS B = A;

   MPO O = raise(L,d);
   MPO O2 = gemm(O,O);

   MPS OA = gemv(O,A);

   clean(OA);

   A = gemv(O,OA);

   clean(A);
   
   MPS OB = gemv(O2,B);

   clean(OB);

   return 0;

}
