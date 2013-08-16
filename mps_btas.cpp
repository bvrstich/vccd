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
using namespace mps;

int main(void){

   cout.precision(10);
   srand(time(NULL));

   //lenght of the chain
   int L = 20;

   //number of particles
   int n_u = 10;
   int n_d = 10;

   Qshapes<Quantum> qp;
   physical(qp);

   MPS A = random(L,Quantum(n_u,n_d),qp,20);
   compress<3>(A,mps::Left,100);
   clean(A);
   normalize(A);

   MPS B = random(L,Quantum(n_u,n_d),qp,20);
   compress<3>(B,mps::Left,100);
   clean(B);
   normalize(B);

   DArray<2> t(10,10);
   t.generate(rgen);

   MPO O = T1(t);

   cout << inprod(mps::Left,A,O,B) << "\t" << inprod(mps::Right,A,O,B) << endl;

   return 0;

}
