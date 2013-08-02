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

   MPO O = ising(L,d,1.0,1.0);

   print(O);
/*
   MPS A = create(L,d,Quantum::zero(),10);

   MPS OA = gemv(O,A);

   print(OA);
*/
   return 0;

}
