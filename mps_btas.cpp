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

   Quantum qt(4);

   MPS A = create(L,qt,d);
   MPS B = create(L,qt,d);

   MPO O = ising(L,d,-1.0,1.0);

   cout << dot(A,B) << "\t" << inprod(A,O,B) << endl;

   return 0;

}
