#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstdlib>

using std::cout;
using std::endl;
using std::ofstream;
using std::ifstream;

#include "SpinQuantum.h"
namespace btas { typedef SpinQuantum Quantum; }; // Defined as default quantum number class

#include "include.h"

using namespace btas;

int main(void){

   cout.precision(10);
   srand(time(NULL));

   //lenght of the chain
   int L = 4;

   //physical dimension
   int d = 2;

   SDArray<3> tmp;

   TVector<Dshapes,3> dshape;

   dshape[0].push_back(1);
   dshape[0].push_back(2);
   dshape[0].push_back(1);

   dshape[1].push_back(1);
   dshape[1].push_back(1);
   dshape[1].push_back(1);

   dshape[2].push_back(2);
   dshape[2].push_back(3);

   tmp.resize(dshape);

   tmp.generate(rgen);

   cout << tmp << endl;

   SDArray<3> tmp2;

   tmp.remove_index(tmp2,0,0);

   cout << tmp2 << endl;

/*
   //MPO O = Sz(L,d);
   MPO O = ising(L,d,1.0,1.0);

   MPS A = create(L,d,Quantum::zero(),10);

   cout << A[0] << endl;
   cout << O[0] << endl;

   MPS OA = gemv(O,A);

   cout << OA[0] << endl;
*/
   return 0;

}
