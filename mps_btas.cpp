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
   int L = 4;

   //physical dimension
   int d = 2;

   //number of particles
   int N = 2;

   MPS A = create(L,d,Quantum(N),10);

   MPO sum = n_loc(L,d,0); 

   for(int i = 1;i < L;++i){

      MPO tmp = add<4>(n_loc(L,d,i),sum);

      sum = tmp;

   }

   print(sum);
   
   cout << inprod(A,sum,A)/nrm2(A) << endl;

   double ward = 0.0;

   for(int i = 0;i < L;++i){

      MPO O = n_loc(L,d,i);

      ward += inprod(A,O,A)/nrm2(A);

   }

   cout << ward << endl;

   return 0;

}
