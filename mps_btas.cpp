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

   cout << endl;
   cout << "A" << endl;
   cout << endl;
   cout << A[0] << endl;

   cout << endl;
   cout << "B" << endl;
   cout << endl;
   cout << B[0] << endl;
   cout << endl;

   QSDArray<3> C;

   QSDjoin(A[0],B[0],C);

   cout << endl;
   cout << "C" << endl;
   cout << endl;
   cout << C << endl;

   return 0;

}
