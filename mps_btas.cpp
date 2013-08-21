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
   int L = 6;

   //number of particles
   int n_u = 3;
   int n_d = 3;

   Qshapes<Quantum> qp;
   physical(qp);

   std::vector<int> occ(L);

   for(int i = 0;i < n_u;++i)
      occ[i] = 3;

   for(int i = n_u;i < L;++i)
      occ[i] = 0;

   MPS<Quantum> A = product_state(L,qp,occ);
   MPS<Quantum> B = create(L,Quantum(n_u,n_d),qp,20,rgen);
   compress(B,mps::Left,100);
   normalize(B);

   DArray<4> t(n_u,n_u,n_u,n_u);
   t.generate(rgen);

   double norm = sqrt(Ddot(t,t));
   Dscal(1.0/norm,t);

   cout << Ddot(t,t) << endl;

   MPO<Quantum> O = T2<Quantum>(t);

   MPS<Quantum> OA = gemv(O,A);
   cout << dot(mps::Left,OA,OA) << endl;

   compress(O,mps::Right,0);
   compress(O,mps::Left,0);
   
   OA = gemv(O,A);
   cout << dot(mps::Left,OA,OA) << endl;

   compress(OA,mps::Right,0);
   compress(OA,mps::Left,0);

   cout << dot(mps::Left,OA,OA) << endl;

   return 0;

}
