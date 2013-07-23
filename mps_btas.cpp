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
  
   btas::Quantum qt(0);

   MPS mps(L,qt,D);

   return 0;

}
