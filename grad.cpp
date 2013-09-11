#include <iostream>
#include <fstream>
#include <cmath>

using std::endl;
using std::ostream;
using std::ofstream;
using std::ifstream;

#include "include.h"

using namespace btas;
using namespace mps;

namespace grad {

   /**
    * construct an MPO, i.e. a QSDArray<4> which contains a fully contracted <B|O|A> with the O[i] missing.
    * @param rol MPS of length L - 1 containing the left renormalized blocks
    * @param ror MPS of length L - 1 containing the right renormalized blocks
    * @param A input MPS
    * @param B input MPS
    */
   MPO<Quantum> construct(const MPS<Quantum> &rol, const MPS<Quantum> &ror,const MPS<Quantum> &A,const MPS<Quantum> &B){

      MPO<Quantum> grad(A.size());

      return grad;

   }

   /**
    * check if the renormalized operators are correctly constructed
    */
   void check(const MPO<Quantum> &grad,const MPO<Quantum> &O){

   }

}
