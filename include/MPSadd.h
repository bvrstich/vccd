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

#include "btas/QSPARSE/QSDArray.h"

namespace btas{

   void QSDjoin(const QSDArray<3> &,const QSDArray<3> &,QSDArray<3> &);

   void QSDjoin_ledge(const QSDArray<3> &,const QSDArray<3> &,QSDArray<3> &);

   void QSDjoin_redge(const QSDArray<3> &,const QSDArray<3> &,QSDArray<3> &);

}

#endif 
