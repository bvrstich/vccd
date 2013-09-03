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
   int L = 8;

   //number of particles
   int n_u = 3;
   int n_d = 3;

   int no = n_u;
   int nv = L - no;

   Ostate::construct_oplist(L);

   Qshapes<Quantum> qp;
   physical(qp);

   MPS<Quantum> A =  mps::create(L,Quantum(n_u,n_d),qp,40,rgen); 
   mps::compress(A,mps::Right,100);
   mps::compress(A,mps::Left,100);
   mps::normalize(A);

   MPS<Quantum> B =  mps::create(L,Quantum(n_u,n_d),qp,40,rgen); 
   mps::compress(B,mps::Right,100);
   mps::compress(B,mps::Left,100);
   mps::normalize(B);

   DArray<4> t(no,no,nv,nv);
   t.generate(rgen);

   for(int i = 0;i < no;++i)
      for(int j = 0;j < no;++j)
         for(int a = 0;a < nv;++a)
            for(int b = 0;b < nv;++b)
               t(i,j,a,b) = t(j,i,b,a);


   MPO<Quantum> T2op = T2<Quantum>(t);
   mps::compress(T2op,mps::Right,0);
   mps::compress(T2op,mps::Left,0);

   MPO<Quantum> T2op_test = T2_test<Quantum>(t);
   mps::compress(T2op_test,mps::Right,0);
   mps::compress(T2op_test,mps::Left,0);

   cout << inprod(mps::Left,A,T2op,B) << endl;
   cout << inprod(mps::Left,B,T2op,A) << endl;
   cout << inprod(mps::Left,A,T2op_test,B)/2.0 << endl;
   cout << inprod(mps::Left,B,T2op_test,A)/2.0 << endl;

   return 0;

}
