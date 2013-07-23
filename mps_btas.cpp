#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

using std::cout;
using std::endl;
using std::ofstream;
using std::ifstream;
using std::vector;

#include "SpinQuantum.h"
namespace btas { typedef SpinQuantum Quantum; };

#include "btas/QSDArray.h"

#include "MPS.h"

int main(void){

   cout.precision(10);

   //lenght of the chain
   int L = 10;

   //physical dimension
   int d = 2;

   //max virtual dimension
   int D = 10;
  
   MPS mps(L);

   btas::Quantum qt(4);

   mps.initialize(qt,D);

   return 0;

}
