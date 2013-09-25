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
   std::vector<int> cutoff(no);
   
   for(int i = 0;i < no;++i)
      cutoff[i] = 0;

   //read in the mp2 guess
   DArray<4> t;
   std::ifstream fin("input/Be/cc-pVDZ/mp2.in");
   boost::archive::binary_iarchive iar(fin);
   iar >> t;
*/
   DArray<4> t(no,no,nv,nv);

   for(int i = 0;i < no;++i)
      for(int j = 0;j < no;++j)
         for(int a = 0;a < nv;++a)
            for(int b = 0;b < nv;++b){

               double value = rgen();

               t(i,j,a,b) = value;
               t(j,i,b,a) = value;

            }

   MPO<Quantum> T = T2<Quantum>(t,true);
   compress(T,mpsxx::Left,0);
   compress(T,mpsxx::Right,0);

   MPS<Quantum> A = create(L,Quantum(n_u,n_d),qp,20,rgen);
   compress(A,mpsxx::Left,100);
   MPS<Quantum> B = create(L,Quantum(n_u,n_d),qp,20,rgen);
   compress(B,mpsxx::Left,100);

   cout << inprod(mpsxx::Left,A,T,B) << endl;

   MPS<Quantum> rol = ro::construct(mpsxx::Left,A,T,B);
   MPS<Quantum> ror = ro::construct(mpsxx::Right,A,T,B);

   MPO<Quantum> grad = grad::construct(rol,ror,A,B);

   T2_2_mpo list(no,nv);

   cout << list << endl;
/*
   int i = 2;
   int a = 1;

   cout << list.get(grad,i,a) << endl;

   //E^0_0
   MPO<Quantum> E_op = E<Quantum>(L,a+no,i,1.0);

   cout << inprod(mpsxx::Left,A,E_op,B) << endl;
*/
   return 0;

}
