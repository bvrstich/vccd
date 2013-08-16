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

   for(int i = 0;i < n_u;++i)//double occupied beneath n_u
      occ[i] = 3;

   for(int i = n_u;i < L;++i)//empty after
      occ[i] = 0;

   MPS<Quantum> A = create(L,Quantum(n_u,n_d),qp,20,rgen);
   
   compress(A,mps::Left,100);
   clean(A);
   normalize(A);

   DArray<2> t(10,10);
   t.generate(rgen);

   double norm = sqrt(2.0)*Dnrm2(t);
   Dscal(1.0/norm,t);

   MPO<Quantum> O = T1<Quantum>(t);
   MPS<Quantum> OA = gemv(O,A);

   cout << dot(mps::Left,OA,OA) << endl;

   return 0;

}
