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
   int L = 14;

   //number of particles
   int n_u = 2;
   int n_d = 2;

   int no = n_u;
   int nv = L - no;

   Ostate::construct_oplist(L);

   Qshapes<Quantum> qp;
   physical(qp);

   std::vector<int> order(L);

   ifstream ord_in("input/Be/cc-pVDZ/order.in");

   for(int i = 0;i < L;++i)
      ord_in >> i >> order[i];

   std::vector<double> e(L);

   ifstream ener_in("input/Be/cc-pVDZ/ener.in");

   for(int i = 0;i < L;++i)
      ener_in >> i >> e[i];

   DArray<2> t(L,L);
   read_oei("input/Be/cc-pVDZ/OEI.in",t,order);
 
   DArray<4> V(L,L,L,L);
   read_tei("input/Be/cc-pVDZ/TEI.in",V,order);

   MPO<Quantum> qc = qcham<Quantum>(t,V);
   compress(qc,mps::Right,0);
   compress(qc,mps::Left,0);

   std::vector<int> occ(L);

   for(int i = 0;i < no;++i)
      occ[i] = 3;

   for(int i = no;i < L;++i)
      occ[i] = 0;

   MPS<Quantum> hf = product_state(L,qp,occ);

   cout << endl;
   cout << "Hartree-Fock " << inprod(mps::Left,hf,qc,hf) << endl;
   cout << endl;

   DArray<4> t2(no,no,nv,nv);
   fill_mp2(t2,V,e);

   vccd::steepest_descent(t2,qc,hf);

   return 0;

}
