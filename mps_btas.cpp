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

   //number of particles
   int n_u = 2;
   int n_d = 2;

   MPS A = create(L,Quantum(n_u,n_d),100);
   compress<3>(A,true,0,true);
   clean(A);

   MPS B = create(L,Quantum(n_u,n_d),100);
   compress<3>(B,true,0,true);
   clean(B);

   MPO O_cr = creator(L,1,1);//spin down on site 1
   MPO O_an = annihilator(L,1,1);//spin down on site 1

   MPO O = gemm(O_cr,O_an);

   MPO O_cr_2 = creator(L,1,0);//spin down on site 1
   MPO O_an_2 = annihilator(L,1,0);//spin down on site 1

   MPO O_2 = gemm(O_cr_2,O_an_2);

   MPO sum = add<4>(O,O_2);

   MPO N = n_loc(L,1);

   cout << inprod(A,sum,A) << "\t" << inprod(A,N,A) << endl;

   return 0;

}
