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

   Qshapes<Quantum> qp;
   physical(qp);

   std::vector<int> occ(L);

   for(int i = 0;i < n_u;++i)
      occ[i] = 3;

   for(int i = n_u;i < L;++i)
      occ[i] = 0;

   MPS<Quantum> A = product_state(L,qp,occ);

   //i j a b
   DArray<4> t(n_u,n_u,n_u,n_u);

   t = 0.0;

   for(int i = 0;i < n_u;++i)
      for(int a = 0;a < n_u;++a)
         t(i,i,a,a) = rgen();

   for(int i = 0;i < n_u;++i)
      for(int a = 0;a < n_u;++a)
         for(int b = a + 1;b < n_u;++b){

            t(i,i,a,b) = rgen();
            t(i,i,b,a) = t(i,i,a,b);

         }

   for(int i = 0;i < n_u;++i)
      for(int j = i + 1;j < n_u;++j)
         for(int a = 0;a < n_u;++a){

            t(i,j,a,a) = rgen();
            t(j,i,a,a) = t(i,j,a,a);

         }

   for(int i = 0;i < n_u;++i)
      for(int j = i + 1;j < n_u;++j)
         for(int a = 0;a < n_u;++a){

            for(int b = 0;b < a;++b){

               t(i,j,a,b) = rgen();
               t(j,i,a,b) = t(i,j,a,b);

            }

            for(int b = a + 1;b < n_u;++b){

               t(i,j,a,b) = rgen();
               t(j,i,a,b) = t(i,j,a,b);

            }

         }

   double norm = Ddot(t,t);
   cout << norm << endl;

   double tmp = 0.0;

   for(int i = 0;i < n_u;++i)
      for(int a = 0;a < n_u;++a)
         tmp += t(i,i,a,a) * t(i,i,a,a);

   for(int i = 0;i < n_u;++i)
      for(int a = 0;a < n_u;++a)
         for(int b = a + 1;b < n_u;++b)
            tmp += 2 * t(i,i,a,b) * t(i,i,a,b);

   for(int i = 0;i < n_u;++i)
      for(int j = i + 1;j < n_u;++j)
         for(int a = 0;a < n_u;++a)
            tmp += 2 * t(i,j,a,a) * t(i,j,a,a);

   for(int i = 0;i < n_u;++i)
      for(int j = i + 1;j < n_u;++j)
         for(int a = 0;a < n_u;++a)
            for(int b = a + 1;b < n_u;++b)
               tmp += 2.0 * (t(i,j,a,b) - t(i,j,b,a)) * (t(i,j,a,b) - t(i,j,b,a));

   //down up coming in: create up
   for(int i = 0;i < n_u;++i)
      for(int j = i + 1;j < n_u;++j)
         for(int a = 0;a < n_u;++a)
            for(int b = a + 1;b < n_u;++b)
               tmp += t(i,j,a,b) * t(i,j,a,b);

   //down up coming in: create down
   for(int i = 0;i < n_u;++i)
      for(int j = i + 1;j < n_u;++j)
         for(int a = 0;a < n_u;++a)
            for(int b = a + 1;b < n_u;++b)
               tmp += t(i,j,b,a) * t(i,j,b,a);

   //up down coming in: create down
   for(int i = 0;i < n_u;++i)
      for(int j = i + 1;j < n_u;++j)
         for(int a = 0;a < n_u;++a)
            for(int b = a + 1;b < n_u;++b)
               tmp += t(i,j,a,b) * t(i,j,a,b);

   //up down coming in: create up
   for(int i = 0;i < n_u;++i)
      for(int j = i + 1;j < n_u;++j)
         for(int a = 0;a < n_u;++a)
            for(int b = a + 1;b < n_u;++b)
               tmp += t(i,j,b,a) * t(i,j,b,a);

   cout << tmp << endl;

   //Dscal(1.0/norm,t);

   MPO<Quantum> O = T2<Quantum>(t);

   MPS<Quantum> OA = gemv(O,A);
   compress(O,mps::Right,0);
   compress(O,mps::Left,0);

   OA = gemv(O,A);
   cout << dot(mps::Left,OA,OA) << endl;

   compress(OA,mps::Right,0);
   compress(OA,mps::Left,0);

   cout << dot(mps::Left,OA,OA) << endl;

   return 0;

}
