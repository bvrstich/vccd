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

   MPO O = XY(L,d,1.0,1.0);

   MPS A = create(L,d,Quantum::zero(),10);
   MPS B = create(L,d,Quantum::zero(),10);

   compress<3>(A,true,0,true);

   return 0;

}
