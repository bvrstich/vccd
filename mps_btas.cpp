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

   MPS A = create(L,Quantum(2),20);
   MPS B = create(L,Quantum(2),20);

   compress(A,true,50);
   compress(B,true,50);

   MPS AB = add(A,B);

   print(AB);

   return 0;

}
