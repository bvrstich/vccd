#ifndef HAMOP_H
#define HAMOP_H

#include <iostream>
#include <cstdlib>
#include <vector>

using std::ostream;

/**
 * @author Brecht Verstichel
 * @date 24-08-2013
 * This class just collects the information needed for incoming and outgoing states of the qc MPO.
 */
class HamOp {

   public:

      //constructor
      HamOp();

      //copy constructor
      HamOp(const HamOp &);

      //destructor
      virtual ~HamOp();

      static void init(int);

      static void print_states();

      static int gL();

      //!list containing the outgoing quantum numbers for every site for the Quantum chemical Hamiltonian
      static std::vector< Qshapes<Quantum> > qo;

      //!list containing the outgoing operators for every site for the Quantum chemical Hamiltonian
      static std::vector< std::vector< Ostate > > ostates;

   private:

      static int L;
      
};

#endif
