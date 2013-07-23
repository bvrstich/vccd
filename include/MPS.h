#ifndef MPS_H
#define MPS_H

#include <iostream>
#include <iomanip>

using std::cout;
using std::endl;

#include "SpinQuantum.h"
namespace btas { typedef SpinQuantum Quantum; };

#include "btas/QSDArray.h"

/**
 * class which contains an array of L btas::QSDArray's. So an MPS of length L
 */
class MPS
{

   public:

      MPS(int);

      MPS(const MPS &);

      virtual ~MPS();

      void initialize(const btas::Quantum &,int);

      int gL() const;

      const btas::QSDArray<3> &operator[](int i);

   private:

      //!length of the chain
      int L;

      //!array containing the mps's
      btas::QSDArray<3> **mps;

};

#endif
