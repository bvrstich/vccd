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

   MPO O = heisenberg(L,d,1.0,1.0,1.0);
   compress<4>(O,true,0,false);
   compress<4>(O,false,0,false);
   clean(O);

   MPO O1 = ising(L,d,1.0,1.0);
   MPO O2 = XY(L,d,1.0,0.0);

   MPO O12 = add<4>(O1,O2);
   clean(O12);
   compress<4>(O12,true,0,false);
   compress<4>(O12,false,0,false);
   clean(O12);

   MPS A = create(L,d,Quantum::zero(),20);
   compress<3>(A,true,100,true);
   clean(A);

   MPS B = create(L,d,Quantum::zero(),20);
   compress<3>(B,true,100,true);
   clean(B);

   MPS OA = gemv(O,A);
   MPS O12A = gemv(O12,A);

   cout << dot(OA,B) << "\t" << dot(O12A,B) << endl;

   return 0;

}
