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
   int L = 20;

   //physical dimension
   int d = 2;

   MPO O1 = ising(L,d,-1.0,1.0);
   MPO O2 = ising(L,d,-1.0,1.0);

   MPO O3 = gemm(O1,O2); 

   return 0;

}
