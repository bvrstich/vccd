#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

using std::cout;
using std::endl;
using std::ofstream;
using std::ifstream;

#include "include.h"

using namespace btas;

int main(void){

   cout.precision(5);
   srand(time(NULL));

   //lenght of the chain
   int L = 20;

   //physical dimension
   int d = 2;

   MPS A = create(L,Quantum(0),20);

   compress(A,true,20);

   for(int i = 1;i < L;++i){

      if(A[i-1].dshape(2) != A[i].dshape(0)){
         cout << A[i-1].qshape() << endl;
         cout << A[i-1].dshape() << endl;
         cout << A[i].qshape() << endl;
         cout << A[i].dshape() << endl;
      }
   }

   return 0;

}
