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

   int no = n_u;
   int nv = L - no;

   Ostate::construct_oplist(L);
   Ostate::print_oplist();

   Qshapes<Quantum> qp;
   physical(qp);

   //i j a b
   DArray<2> t(L,L);
   t.generate(rgen);

   MPO<Quantum> O = one_body<Quantum>(t);

   return 0;

}
