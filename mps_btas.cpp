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
   int L = 100;

   //physical dimension
   int d = 3;

   MPS A = create(L,d,Quantum::zero(),10);
   compress(A,true,100);

   clean(A);

   compress(A,false,0);

   for(int i = 0;i < L;++i){

      cout << A[i].qshape() << endl;
      cout << A[i].dshape() << endl;

   }


      

   return 0;

}
