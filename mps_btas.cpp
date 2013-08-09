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

#include "include.h"

using namespace btas;

int main(void){

   cout.precision(10);
   srand(time(NULL));

   //lenght of the chain
   int L = 20;

   //physical dimension
   int d = 2;

   //number of particles
   int N = 11;

   MPS A = create(L,d,Quantum(N),20);
   compress<3>(A,true,100,true);
   clean(A);

   MPS B = create(L,d,Quantum(N),20);
   compress<3>(B,true,100,true);
   clean(B);

   MPO O = hopping(L,d);

   MPS OA = gemv(O,A);
   compress<3>(OA,true,0,false);
   clean(OA);

   cout << dot(OA,B) << endl;

   MPO O_cr = creator(L,d,0);
   MPO O_an = annihilator(L,d,1);
   MPO sum1 = gemm(O_cr,O_an);
   MPO tmp,tmp2;

   for(int i = 1;i < L - 1;++i){

      O_cr = creator(L,d,i);
      O_an = annihilator(L,d,i + 1);

      tmp = gemm(O_cr,O_an);

      tmp2 = add<4>(tmp,sum1);

      sum1 = tmp2;

   }

   O_cr = creator(L,d,1);
   O_an = annihilator(L,d,0);
   MPO sum2 = gemm(O_cr,O_an);

   for(int i = 1;i < L - 1;++i){

      O_cr = creator(L,d,i + 1);
      O_an = annihilator(L,d,i);

      tmp = gemm(O_cr,O_an);

      tmp2 = add<4>(tmp,sum2);

      sum2 = tmp2;

   }
   
   MPO sum = add<4>(sum1,sum2);
   clean(sum);
   compress<4>(sum,true,100,false);
   clean(sum);
   compress<4>(sum,false,100,false);
   clean(sum);

   cout << inprod(A,sum,B) << endl;

   for(int i = 0;i < L;++i){

      cout << endl;
      cout << "site " << i << endl;
      cout << endl;
      cout << O[i].qshape() << endl;
      cout << O[i].dshape() << endl;
      cout << sum[i].qshape() << endl;
      cout << sum[i].dshape() << endl;
      cout << endl;

   }
   return 0;

}
