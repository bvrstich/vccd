#ifndef LS_H
#define LS_H

#include <iostream>
#include <cstdlib>
#include <vector>

using std::ostream;

using namespace mpsxx;
using namespace btas;

class e_eMPS;

/**
 * @author Brecht Verstichel
 * @date 11-09-2013
 * this namespace contains some functions which create the objects from which ls can be calculated efficiently
 */
namespace ls {

      //constructor
      template<class Q>
         void construct(DArray<2> &,DArray<2> &,const MPO<Q> &qc,const MPS<Q> &hf,const eMPS &,const e_eMPS &);

      double eval(const DArray<2> &,const DArray<2> &,double alpha);

}

#endif
