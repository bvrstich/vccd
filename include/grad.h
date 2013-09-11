#ifndef RENOP_H
#define RENOP_H

#include <iostream>
#include <cstdlib>
#include <vector>

using std::ostream;

using namespace mps;
using namespace btas;

/**
 * @author Brecht Verstichel
 * @date 11-09-2013
 * this namespace contains some functions which create the objects from which gradients can be calculated efficiently
 */
namespace grad {

      //constructor
      MPO<Quantum> construct(const MPS<Quantum> &rol,const MPS<Quantum> &ror,const MPS<Quantum> &A,const MPS<Quantum> &B);

      void check(const MPO<Quantum> &grad,const MPO<Quantum> &O);

}

#endif
