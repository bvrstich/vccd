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

int main(void){

   cout.precision(10);
   srand(time(NULL));

   //lenght of the chain
   int L = 8;

   //number of particles
   int n_u = 4;
   int n_d = 4;

   int no = n_u;
   int nv = L - no;

   Ostate::construct_oplist(L);

   Qshapes<Quantum> qp;
   physical(qp);
/*
   //make the HF state
   std::vector<int> occ(L);

   for(int i = 0;i < no;++i)
      occ[i] = 3;

   for(int i = no;i < L;++i)
      occ[i] = 0;

   MPS<Quantum> hf = product_state(L,qp,occ);

   //load the qc hamiltonian
   MPO<Quantum> qc(L);
   load_mpx(qc,"input/Be/cc-pVDZ/MPO/qcham");

   //hartree fock energy
   cout << inprod(mpsxx::Left,hf,qc,hf) << endl;

   //the cutoff vector for the exponential
   std::vector<int> cutoff(1);

   cutoff[0] = 20;
   //cutoff[1] = 10;

   //read in the mp2 guess
   DArray<4> t(no,no,nv,nv);

   std::ifstream fin("input/Be/cc-pVDZ/mp2.in");
   boost::archive::binary_iarchive iar(fin);
   iar >> t;

   vccd::conjugate_gradient(t,qc,hf,cutoff);
   //vccd::steepest_descent(t,qc,hf,cutoff);
*/
   DArray<2> t(L,L);

   for(int i = 0;i < L;++i)
      for(int j = i;j < L;++j){

         double value = rgen();

         t(i,j) = value;
         t(j,i) = value;

      }

   DArray<4> V(L,L,L,L);

   for(int i = 0;i < L;++i)
      for(int j = 0;j < L;++j)
         for(int k = 0;k < L;++k)
            for(int l = 0;l < L;++l){

               double value = rgen();

               V(i,j,k,l) = value;
               V(j,i,l,k) = value;
               V(k,j,i,l) = value;
               V(j,k,l,i) = value;
               V(i,l,k,j) = value;
               V(l,i,j,k) = value;
               V(k,l,i,j) = value;
               V(l,k,j,i) = value;

            }


   MPS<Quantum> A = create(L,Quantum(n_u,n_d),qp,20,rgen());

   compress(A,mpsxx::Left,0);
   compress(A,mpsxx::Right,0);

   normalize(A);

   MPS<Quantum> B = create(L,Quantum(n_u,n_d),qp,20,rgen());

   compress(B,mpsxx::Left,0);
   compress(B,mpsxx::Right,0);

   normalize(B);

   MPO<Quantum> qc = qcham<Quantum>(t,V);

   compress(qc,mpsxx::Left,0);
   compress(qc,mpsxx::Right,0);

   cout << inprod(mpsxx::Left,A,qc,B) << endl;
   cout << inprod(mpsxx::Left,B,qc,A) << endl;

   MPO<Quantum> qc_new = qcham_new<Quantum>(t,V,false);

   compress(qc_new,mpsxx::Left,0);
   compress(qc_new,mpsxx::Right,0);

   for(int i = 0;i < L;++i){

      cout << endl;
      cout << "site " << i << endl;
      cout << endl;
      cout << qc_new[i].qshape() << endl;
      cout << qc_new[i].dshape() << endl;
      cout << endl;

   }

   return 0;

}
