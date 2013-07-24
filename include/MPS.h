#ifndef MPS_H
#define MPS_H

#include <iostream>
#include <iomanip>

using std::cout;
using std::endl;
using std::ostream;

#include "Tensor.h"

using namespace btas;

/**
 * class which contains an array of L Tensor objects. So an MPS of length L. 
 * In this class all of the MPS related things are calculated, such as the quantumnumbers and dimensions of the blocks
 * present in the Tensors, canonicalization of the chain, normalization etc..
 */
class MPS
{
   /**
    * Output stream operator overloaded
    * @param output The stream to which you are writing (e.g. cout)
    * @param mps_p the mps you want to print
    */
   friend ostream &operator<<(ostream &output,MPS &mps_p);

   public:

      MPS(int,const Quantum &,int);

      MPS(const MPS &);

      virtual ~MPS();

      void initialize();

      int gL() const;

      const Tensor &operator[](int i) const;

      const Quantum &gqt() const;

      int gD() const;

      void canonicalize(bool);

   private:

      //!length of the chain
      int L;

      //!maximal dimension of the symmetryblocks
      int D;

      //!total quantumnumber of the chain
      Quantum *qt;

      //!array containing the mps's
      Tensor **mps;

};

#endif
