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

   MPS A = create(L,Quantum(2),10);
   MPS B = create(L,Quantum(2),20);

   MPS AB = axpy(1.0,A,B);

   print(AB);

   return 0;

}
