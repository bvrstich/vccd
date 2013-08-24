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
   int L = 6;

   //number of particles
   int n_u = 3;
   int n_d = 3;

   int no = n_u;
   int nv = L - no;

   Qshapes<Quantum> qp;
   physical(qp);

   //i j a b
   DArray<4> t(no,no,nv,nv);
   t.generate(rgen);

   t = 0.0;

   for(int i = 0;i < no;++i)
      for(int a = 0;a < nv;++a)
         t(i,i,a,a) = rgen();

   for(int i = 0;i < no;++i)
      for(int a = 0;a < nv;++a)
         for(int b = a + 1;b < nv;++b){

            t(i,i,a,b) = rgen();
            t(i,i,b,a) = t(i,i,a,b);

         }

   for(int i = 0;i < no;++i)
      for(int j = i + 1;j < no;++j)
         for(int a = 0;a < nv;++a){

            t(i,j,a,a) = rgen();
            t(j,i,a,a) = t(i,j,a,a);

         }

   for(int i = 0;i < no;++i)
      for(int j = i + 1;j < no;++j)
         for(int a = 0;a < nv;++a){

            for(int b = 0;b < a;++b){

               t(i,j,a,b) = rgen();
               t(j,i,a,b) = t(i,j,a,b);

            }

            for(int b = a + 1;b < nv;++b){

               t(i,j,a,b) = rgen();
               t(j,i,a,b) = t(i,j,a,b);

            }

         }

   double norm = Ddot(t,t);
   cout << norm << endl;

   double tmp = 0.0;

   for(int i = 0;i < no;++i)
      for(int a = 0;a < nv;++a)
         tmp += t(i,i,a,a) * t(i,i,a,a);

   for(int i = 0;i < no;++i)
      for(int a = 0;a < nv;++a)
         for(int b = a + 1;b < nv;++b)
            tmp += 2 * t(i,i,a,b) * t(i,i,a,b);

   for(int i = 0;i < no;++i)
      for(int j = i + 1;j < no;++j)
         for(int a = 0;a < nv;++a)
            tmp += 2 * t(i,j,a,a) * t(i,j,a,a);

   for(int i = 0;i < no;++i)
      for(int j = i + 1;j < no;++j)
         for(int a = 0;a < nv;++a)
            for(int b = a + 1;b < nv;++b)
               tmp += 2.0 * (t(i,j,a,b) - t(i,j,b,a)) * (t(i,j,a,b) - t(i,j,b,a));

   for(int i = 0;i < no;++i)
      for(int j = i + 1;j < no;++j)
         for(int a = 0;a < nv;++a)
            for(int b = a + 1;b < nv;++b)
               tmp += 2*t(i,j,a,b) * t(i,j,a,b) + 2*t(i,j,b,a) * t(i,j,b,a);

   cout << tmp << endl;

   //Dscal(1.0/norm,t);

   MPO<Quantum> O = T2<Quantum>(t);
   compress(O,mps::Right,0);
   compress(O,mps::Left,0);

   std::vector<int> occ(L);

   for(int i = 0;i < no;++i)
      occ[i] = 3;

   for(int i = no;i < L;++i)
      occ[i] = 0;

   MPS<Quantum> A = product_state(L,qp,occ);

   MPS<Quantum> OA = gemv(O,A);
   cout << dot(mps::Left,OA,OA) << endl;

   return 0;

}
