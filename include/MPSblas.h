#ifndef _BTAS_MPSBLAS_H
#define _BTAS_MPSBLAS_H 1

#include <iostream>
#include <iomanip>

using std::cout;
using std::endl;
using std::ostream;

#include "SpinQuantum.h"
namespace btas { typedef SpinQuantum Quantum; };

#include "btas/QSDArray.h"

namespace btas{

   //!typedefine MPS as an std::vector< QSDArray<3> > 
   typedef std::vector< QSDArray<3> > MPS;

   //some function definitions on MPS's
   MPS create(int,const Quantum &qt,int);

   double rgen();

   void print(const MPS &);

   void copy(const MPS &,MPS &);

   void scal(double,MPS &);

   void compress(MPS &,bool,int);

   void axpy(double,const MPS &,MPS &);

   QSDArray<2> dot(const MPS &,const MPS &);

   double nrm2(const MPS &);

}

#endif 
