#ifndef _BTAS_MPSADD_H
#define _BTAS_MPSADD_H 1

#include <iostream>
#include <iomanip>

using std::cout;
using std::endl;
using std::ostream;

//some functions needed for the addiotion of two MPS's

#include "SpinQuantum.h"
namespace btas { typedef SpinQuantum Quantum; };

#include "btas/QSDArray.h"

namespace btas{

   void qindex(const QSDArray<3> &,int,Qshapes &);

   void Djoin(const DArray<3> &,const DArray<3> &,DArray<3> &);

   void QSDjoin(const QSDArray<3> &,const QSDArray<3> &,QSDArray<3> &);

}

#endif 
