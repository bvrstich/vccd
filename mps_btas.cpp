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
   int n_u = 9;
   int n_d = 9;

   int no = n_u;
   int nv = L - no;

   Ostate::construct_oplist(L);

   Qshapes<Quantum> qp;
   physical(qp);

   //i j a b
   DArray<4> t(no,no,nv,nv);
   t.generate(rgen);

   ifstream in("t.in");

   for(int i = 0;i < no;++i)
      for(int j = 0;j < no;++j)
         for(int a = 0;a < nv;++a)
            for(int b = 0;b < nv;++b)
               in >> i >> j >> a >> b >> t(i,j,a,b);

   MPO<Quantum> O = T2<Quantum>(t);

   for(int i = 0;i < L;++i){
      cout << endl;
      cout << "site " << i << endl;
      cout << endl;
      cout << O[i] << endl;
   }


   return 0;

}
