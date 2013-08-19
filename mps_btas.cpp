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

/**
 * simple random number generator
 */
double rgen() { return 2.0*(static_cast<double>(rand())/RAND_MAX) - 1.0; }

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

   std::vector<int> occ(L);

   MPS<Quantum> A = create(L,Quantum(n_u,n_d),qp,20,rgen);
   compress(A,mps::Left,100);
   MPS<Quantum> B = create(L,Quantum(n_u,n_d),qp,20,rgen);
   compress(B,mps::Left,100);

   MPS<Quantum> AB = add(A,B);

   MPS<Quantum> C = create(L,Quantum(n_u,n_d),qp,20,rgen);

   cout << dot(mps::Left,A,C) + dot(mps::Left,B,C) << "\t" << dot(mps::Left,AB,C) << endl;

   return 0;

}
