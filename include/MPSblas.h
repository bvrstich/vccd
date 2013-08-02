#ifndef _BTAS_MPSBLAS_H
#define _BTAS_MPSBLAS_H 1

#include <iostream>
#include <iomanip>

using std::cout;
using std::endl;
using std::ostream;

#include "SpinQuantum.h"
namespace btas { typedef SpinQuantum Quantum; };

#include "btas/QSPARSE/QSDArray.h"

namespace btas{

   //!typedefine MPS as an std::vector< QSDArray<3> > 
   typedef std::vector< QSDArray<3> > MPS;

   //!typedefine MPO as an std::vector< QSDArray<4> > 
   typedef std::vector< QSDArray<4> > MPO;

   //some function definitions on MPS's
   MPS create(int,const Quantum &qt,int);

   double rgen();

   void print(const MPS &);

   void print(const MPO &);

   void copy(const MPS &,MPS &);

   void scal(double,MPS &);

   void compress(MPS &,bool,int);

   MPS add(const MPS &,const MPS &);

   double dot(const MPS &,const MPS &);

   double nrm2(const MPS &);

   double dist(const MPS &,const MPS &);

   double inprod(const MPS &,const MPO &,const MPS &);

   MPS gemv(const MPO &,const MPS &);

   MPO gemm(const MPO &,const MPO &);

}

#endif 
