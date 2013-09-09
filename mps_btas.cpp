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

   MPS<Quantum> X = create(L,Quantum(n_u,n_d),qp,100,rgen);
   normalize(X);
   MPS<Quantum> Y = create(L,Quantum(n_u,n_d),qp,100,rgen);
   normalize(Y);

   DArray<2> tA(L,L);

   for(int i = 0;i < L;++i)
      for(int j = 0;j < L;++j){

         double value = rgen();

         tA(i,j) = value;
         tA(j,i) = value;

      }


   MPO<Quantum> A = one_body<Quantum>(tA);
   compress(A,mps::Right,0);
   compress(A,mps::Left,0);

   for(int i = 0;i < L;++i)
      for(int j = 0;j < L;++j){

         double value = rgen();

         tA(i,j) = value;
         tA(j,i) = value;

      }

   MPO<Quantum> B = one_body<Quantum>(tA);
   compress(B,mps::Right,0);
   compress(B,mps::Left,0);

   for(int i = 0;i < L;++i)
      for(int j = 0;j < L;++j){

         double value = rgen();

         tA(i,j) = value;
         tA(j,i) = value;

      }

   MPO<Quantum> C = one_body<Quantum>(tA);
   compress(C,mps::Right,0);
   compress(C,mps::Left,0);

   MPO<Quantum> C_copy(C);

   double alpha = 0.62734;
   double beta = 0.12349;

   MPO<Quantum> AB = A*B;

   gemm(alpha,A,B,beta,C);

   cout << alpha * inprod(mps::Left,X,AB,Y) + beta * inprod(mps::Left,X,C_copy,Y) << "\t" << inprod(mps::Left,X,C,Y) << endl;

   /*
   //make the HF state
   std::vector<int> occ(L);

   for(int i = 0;i < no;++i)
   occ[i] = 3;

   for(int i = no;i < L;++i)
   occ[i] = 0;

   MPS<Quantum> hf = product_state(L,qp,occ);
   MPS<Quantum> HA = qc*hf;
    */
   return 0;

}
