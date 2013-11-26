#ifndef GRAD_H
#define GRAD_H

#include <iostream>
#include <cstdlib>
#include <vector>

using std::ostream;

using namespace mpsxx;
using namespace btas;

/**
 * @author Brecht Verstichel
 * @date 11-09-2013
 * this namespace contains some functions which create the objects from which gradients can be calculated efficiently
 */
namespace grad {

      //constructor
      MPO<Quantum> construct(int,int,double,const MPS<Quantum> &);

      void get_op(int site,int opnum,const QSDArray<3> &A,QSDArray<5> &);

      void contract_op(int site,int opnum,const QSDArray<2> &op,const QSDArray<5> &E,QSDArray<4> &G);

}

#endif
