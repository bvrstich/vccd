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
using namespace mpsxx;

int main(int argc,char *argv[]){

   char *dirpath = argv[1];

   cout.precision(10);
   srand(time(NULL));

   //lenght of the chain
   int L = 14;

   //number of particles
   int n_u = 5;
   int n_d = 5;

   int no = n_u;
   int nv = L - no;

   Ostate::construct_oplist(L);

   Qshapes<Quantum> qp;
   physical(qp);

   //make the HF state
   std::vector<int> occ(L);

   for(int i = 0;i < no;++i)
      occ[i] = 3;

   for(int i = no;i < L;++i)
      occ[i] = 0;

   //hf energies
   std::vector<double> e;

   char enerfile[100];
   sprintf(enerfile,"%sener.in",dirpath);

   std::ifstream ener_in(enerfile);

   int i;
   double value;

   while(ener_in >> i >> value)
      e.push_back(value);

   MPS<Quantum> hf = product_state(L,qp,occ);

   std::vector<int> order;

   char orderfile[100];
   sprintf(orderfile,"%sorder.in",dirpath);

   std::ifstream ord_in(orderfile);

   int j;

   while(ord_in >> i >> j)
      order.push_back(j);

   char oeifile[100];
   sprintf(oeifile,"%sOEI.in",dirpath);

   //construct the qc hamiltonian
   DArray<2> K(L,L);
   read_oei(oeifile,K,order);

   char teifile[100];
   sprintf(teifile,"%sTEI.in",dirpath);

   //construct the qc hamiltonian
   DArray<4> V(L,L,L,L);
   read_tei(teifile,V,order);

   MPO<Quantum> qc = qcham<Quantum>(K,V,false);

   cout << compress(qc,mpsxx::Left,0) << endl;
   cout << compress(qc,mpsxx::Right,0) << endl;

   print_dim(qc);

   //hartree fock energy
   cout << inprod(mpsxx::Left,hf,qc,hf) << endl;

   //construct the mp2 guess
   DArray<4> t(no,no,nv,nv);

   fill_mp2(t,V,e);
   
   //solve
   //vccd::solve(t,qc,hf,e,0,0);
   vccd::conjugate_gradient(t,qc,hf,e,100,no);

   return 0;

}
