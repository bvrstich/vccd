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

   ifstream in_order("input/Be/cc-pVDZ/order.in");

   for(int i = 0;i < L;++i)
      in_order >> i >> order[i]; 

   DArray<2> t(L,L);
   t = 0;
   read_oei("input/Be/cc-pVDZ/OEI.in",t,order);

   DArray<4> V(L,L,L,L);
   V = 0.0;
   read_tei("input/Be/cc-pVDZ/TEI.in",V,order);

   std::vector<double> e(L);

   ifstream in_ener("input/Be/cc-pVDZ/ener.in");

   for(int i = 0;i < L;++i)
      in_ener >> i >> e[i]; 

   //construct quantum chemistry MPO
   MPO<Quantum> qc = qcham<Quantum>(t,V);
   compress(qc,mps::Right,0);
   compress(qc,mps::Left,0);

   //input HF state
   std::vector<int> occ(L);

   //bra
   for(int i = 0;i < L;++i)
      occ[i] = 0;

   for(int i = 0;i < no;++i)
      occ[i] = 3;

   MPS<Quantum> A = product_state(L,qp,occ);

   /*
      double hf = inprod(mps::Left,A,qc,A);

   //T2 operator: fill with MP2 input
   DArray<4> t2(no,no,nv,nv);
   fill_mp2(t2,V,e);

   //construct T2 MPO
   MPO<Quantum> T2_op = T2<Quantum>(t2);
   compress(T2_op,mps::Right,0);
   compress(T2_op,mps::Left,0);

   MPS<Quantum> eTA = exp(T2_op,A,no);
   cout << inprod(mps::Left,A,qc,eTA) << endl;

   normalize(eTA);
   cout << dot(mps::Left,eTA,eTA) << "\t" << inprod(mps::Left,eTA,qc,eTA) << endl;

   cout << eTA << endl;
    */
   return 0;

}
