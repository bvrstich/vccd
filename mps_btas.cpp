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
   HamOp::init(L);

   //hf energies
   std::vector<double> e;

   char enerfile[100];
   sprintf(enerfile,"%sener.in",dirpath);

   std::ifstream ener_in(enerfile);

   int i;
   double value;

   while(ener_in >> i >> value)
      e.push_back(value);

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
   //read_oei(oeifile,K,order);
   random_oei(K);

   char teifile[100];
   sprintf(teifile,"%sTEI.in",dirpath);

   //construct the qc hamiltonian
   DArray<4> V(L,L,L,L);
   //read_tei(teifile,V,order);
   random_tei(V);

   //construct the complementary operators
   Operator::init(K,V);

   //here start the construction of states and stuff..., the actual program
   Qshapes<Quantum> qp;
   physical(qp);

   //make the HF state
   std::vector<int> occ(L);

   for(int i = 0;i < no;++i)
      occ[i] = 3;

   for(int i = no;i < L;++i)
      occ[i] = 0;

   MPS<Quantum> hf = product_state(L,qp,occ);

   MPO<Quantum> qc = qcham<Quantum>(K,V,false);

   cout << compress(qc,mpsxx::Left,0) << endl;
   cout << compress(qc,mpsxx::Right,0) << endl;

   print_dim(qc);

   //hartree fock energy
   cout << inprod(mpsxx::Left,hf,qc,hf) << endl;
/*
   //construct the mp2 guess
   DArray<4> t(no,no,nv,nv);

   fill_mp2(t,V,e);

   MPO<Quantum> T = T2<Quantum>(t,false);

   cout << compress(T,mpsxx::Left,0) << endl;
   cout << compress(T,mpsxx::Right,0) << endl;

   MPS<Quantum> eTA = exp(T,hf,no,100);
   normalize(eTA);

   print_dim(eTA);

   cout << inprod(mpsxx::Left,eTA,qc,eTA) << endl;
 */ 
   
   MPS<Quantum> A = mpsxx::create(L,Quantum(n_u,n_d),qp,40,rgen);

   cout << compress(A,mpsxx::Left,100) << endl;
   cout << compress(A,mpsxx::Right,100) << endl;

   //MPO<Quantum> tmp = ro::construct(mpsxx::Left,eTA,T,qc,eTA);

   ro::construct(mpsxx::Left,A);

   return 0;

}
