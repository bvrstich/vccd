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

   cout.precision(10);
   srand(time(NULL));

   //lenght of the chain
   int L = 20;

   //physical dimension
   int d = 2;

   MPS mps_X = create(L,Quantum(2),10);
   MPS mps_Y = create(L,Quantum(2),10);

   cout << dist(mps_X,mps_Y) << endl;

   return 0;

}
