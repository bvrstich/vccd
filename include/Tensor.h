#ifndef TENSOR_H
#define TENSOR_H

#include <iostream>
#include <iomanip>

using std::cout;
using std::endl;
using std::ostream;

#include "SpinQuantum.h"
namespace btas { typedef SpinQuantum Quantum; };

#include "btas/QSDArray.h"

using namespace btas;

/**
 * This class enherits from and specializes the block-sparse QSDArray<> template class for a regular 3-index tensor for a 1D MPS.
 * It extends on the class by defining special functions specific to the 3-index tensor from a 1D MPS.
 */
class Tensor : public QSDArray<3>
{
   public:

      Tensor();
      
      Tensor(const blitz::TinyVector<Qshapes,3> &,const blitz::TinyVector<Dshapes,3> &);

      Tensor(const Tensor &);

      virtual ~Tensor();

      using QSDArray<3>::operator=;

      void canonicalize(bool);

      void fill_Random();

   private:

};

#endif
