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

      void construct(const MPS_DIRECTION &dir,const MPS<Quantum> &wccd);

      void check(const MPS<Quantum> &ror,const MPS<Quantum> &rol);

      void check(const MPO<Quantum> &ror,const MPO<Quantum> &rol);

      void print_op(int site,int opnum,const QSDArray<2> &op,const QSDArray<3> &A);

      void print_op(int site,int opnum,const QSDArray<2> &op,const QSDArray<4> &A);

      void get_op(int site,int opnum,const QSDArray<3> &A,QSDArray<4> &);

      void calc_op(const QSDArray<2> &op,const QSDArray<4> &A,QSDArray<2> &);

      void read(int site,int opnum,QSDArray<2> &);
      
      void save(int site,int opnum,const QSDArray<2> &);

}

#endif
