#ifndef RO_H
#define RO_H

#include <iostream>
#include <cstdlib>
#include <vector>

using std::ostream;

using namespace mpsxx;
using namespace btas;

class TeMPS;

/**
 * @author Brecht Verstichel
 * @date 10-09-2013
 * this namespace contains the left and right renormalized operators for a certain MPS/MPO inner product.
 */
namespace ro {

      //constructor
      MPS<Quantum> construct(const MPS_DIRECTION &dir,const MPS<Quantum> &,const MPO<Quantum> &,const MPS<Quantum> &);

      MPO<Quantum> construct(const MPS_DIRECTION &dir,const MPS<Quantum> &,const MPO<Quantum> &,const MPO<Quantum> &,const MPS<Quantum> &);

      void check(const MPS<Quantum> &ror,const MPS<Quantum> &rol);

      void check(const MPO<Quantum> &ror,const MPO<Quantum> &rol);

}

#endif
