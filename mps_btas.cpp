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

   MPS A = create(L,d,Quantum::zero(),20);
   MPS B = create(L,d,Quantum::zero(),20);

   //after compressing: clean!
   compress(A,true,100);
   clean(A);

   compress(B,true,100);
   clean(B);

   MPS AB = add(A,B);

   MPS C = create(L,d,Quantum::zero(),20);

   cout << dot(AB,C) << "\t" << dot(A,C) + dot(B,C) << endl;

   return 0;

}
