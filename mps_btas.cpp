#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstdlib>

using std::cout;
using std::endl;
using std::ofstream;
using std::ifstream;

#include "FermiQuantum.h"
namespace btas { typedef FermiQuantum Quantum; }; // Defined as default quantum number class

#include "include.h"

using namespace btas;

int main(void){

   cout.precision(10);
   srand(time(NULL));

   //lenght of the chain
   int L = 20;

   //number of particles
   int n_u = 10;
   int n_d = 9;

   MPS A = HF(L,Quantum(n_u,n_d));

   for(int i = 0;i < L;++i){

      cout << endl;
      cout << A[i].qshape() << endl;
      cout << A[i].dshape() << endl;
      cout << endl;

   }


   return 0;

}
