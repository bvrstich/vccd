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
using namespace mps;

int main(void){

   cout.precision(10);
   srand(time(NULL));

   //lenght of the chain
   int L = 20;

   //number of particles
   int n_u = 10;
   int n_d = 10;

   Qshapes<Quantum> qp;
   physical(qp);

   MPS A = random(L,Quantum(n_u,n_d),qp,20);
   compress<3>(A,true,20,true);

   MPS B = random(L,Quantum(n_u,n_d),qp,20);
   compress<3>(B,true,20,true);

   MPS AB = add<3>(A,B);

   for(int i = 1;i < L;++i){

      if(AB[i].dshape(0) != AB[i - 1].dshape(2)){

         cout << AB[i].dshape(0) << endl;
         cout << AB[i - 1].dshape(2) << endl;

      }

   }


   MPS C = random(L,Quantum(n_u,n_d),qp,20);
   compress<3>(B,true,20,true);

   cout << dot(AB,C) << endl;

   return 0;

}
