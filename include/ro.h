#ifndef RO_H
#define RO_H

#include <iostream>
#include <cstdlib>
#include <vector>

using std::ostream;

using namespace mpsxx;
using namespace btas;

/**
 * @author Brecht Verstichel
 * @date 10-09-2013
 * this namespace contains the left and right renormalized operators for a certain MPS/MPO inner product.
 */
namespace ro {

      void construct(const MPS_DIRECTION &dir,const MPO<Quantum> &T,const MPS<Quantum> &wccd);

      void print_op(const MPS_DIRECTION &dir,int site,int opnum,const QSDArray<2> &op,const QSDArray<4> &,const QSDArray<3> &A);

      void print_op(const MPS_DIRECTION &dir,int site,int opnum,const QSDArray<2> &op,const QSDArray<5> &A);

      void get_op(const MPS_DIRECTION &dir,int site,int opnum,const QSDArray<4> &,const QSDArray<3> &A,QSDArray<5> &);

      void read(const MPS_DIRECTION &dir,int site,int opnum,QSDArray<3> &);
      
      void save(const MPS_DIRECTION &dir,int site,int opnum,const QSDArray<3> &);

}

#endif
