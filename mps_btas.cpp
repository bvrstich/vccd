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
   int L = 14;

   //number of particles
   int n_u = 2;
   int n_d = 2;

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
   DArray<4> t(no,no,nv,nv);

   std::ifstream fin("input/Be/cc-pVDZ/mp2.in");
   boost::archive::binary_iarchive iar(fin);
   iar >> t;

   vccd::conjugate_gradient(t,qc,hf,cutoff);
/*
   MPO<Quantum> T = T2<Quantum>(t,false);

   MPS<Quantum> eTA = exp(T,hf,cutoff);
   normalize(eTA);

   MPS<Quantum> rol = ro::construct(mpsxx::Left,eTA,T,eTA);
   MPS<Quantum> ror = ro::construct(mpsxx::Right,eTA,T,eTA);

   MPO<Quantum> grad = grad::construct(rol,ror,eTA,eTA);

   T2_2_mpo list(no,nv);

   for(int i = 0;i < no;++i){

      for(int a = 0;a < nv;++a)
         for(int b = 0;b < nv;++b){

            double tmp1 = list.get(grad,i,i,b,a) + list.get(grad,i,i,a,b);

            //E^a_i E^b_j
            MPO<Quantum> E_op = E<Quantum>(L,a+no,b+no,i,i,1.0);
            double tmp2 = inprod(mpsxx::Left,eTA,E_op,eTA);

            cout << endl;
            cout << i << "\t" << i << "\t" << a << "\t" << b << endl;
            cout << endl;
            cout << tmp1 << "\t" << tmp2 << endl;

         }

      for(int j = i + 1;j < no;++j)
         for(int a = 0;a < nv;++a)
            for(int b = 0;b < nv;++b){

               double tmp1 = list.get(grad,i,j,a,b);

               //E^a_i E^b_j
               MPO<Quantum> E_op = E<Quantum>(L,a+no,b+no,i,j,1.0);
               double tmp2 = inprod(mpsxx::Left,eTA,E_op,eTA);

               cout << endl;
               cout << i << "\t" << j << "\t" << a << "\t" << b << endl;
               cout << endl;
               cout << tmp1 << "\t" << tmp2 << endl;

            }

   }
*/
   return 0;

}
