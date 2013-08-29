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

   std::vector<int> order = {0,1,5,12,13,9,10,11,2,6,3,7,4,8};

   DArray<2> t(L,L);
   t = 0;
   read_oei("input/Be/cc-pVDZ/OEI.in",t,order);

   DArray<4> V(L,L,L,L);
   V = 0.0;
   read_tei("input/Be/cc-pVDZ/TEI.in",V,order);

   std::vector<double> e(L);

   ifstream input("input/Be/cc-pVDZ/ener.in");

   for(int i = 0;i < L;++i)
      input >> i >> e[i]; 

   MPO<Quantum> O = qcham<Quantum>(t,V);
   compress(O,mps::Right,0);
   compress(O,mps::Left,0);

   std::vector<int> occ(L);

   //bra
   for(int i = 0;i < L;++i)
      occ[i] = 0;

   occ[0] = 3;
   occ[1] = 3;

   MPS<Quantum> A = product_state(L,qp,occ);

   DArray<4> t2(no,no,nv,nv);
   fill_mp2(t2,V,e);

   MPO<Quantum> T2_op = T2<Quantum>(t2);

   compress(T2_op,mps::Right,0);
   compress(T2_op,mps::Left,0);

   MPS<Quantum> T2A = gemv(T2_op,A);

   double mp2 = 0.0;

   for(int i = 0;i < no;++i)
      for(int j = 0;j < no;++j)
         for(int a = no;a < L;++a)
            for(int b = no;b < L;++b)
               mp2 += 2.0 * V(i,j,a,b) * V(i,j,a,b) / ( e[i] + e[j] - e[a] - e[b] ) - V(i,j,a,b) * V(i,j,b,a) / ( e[i] + e[j] - e[a] - e[b] );

   cout << mp2 << endl;
   cout << inprod(mps::Left,A,O,A) + inprod(mps::Left,A,O,T2A) << endl;

   return 0;

}
