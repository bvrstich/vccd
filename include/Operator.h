#ifndef OPERATOR_H
#define OPERATOR_H

#include <iostream>
#include <cstdlib>
#include <vector>

using std::ostream;

/**
 * @author Brecht Verstichel
 * @date 24-08-2013
 * This class contains the different operators needed in the construction of the renormalized operators for the gradient
 */
class Operator {

   public:

      //constructor
      Operator();

      //copy constructor
      Operator(const Operator &);

      //destructor
      virtual ~Operator();

      static void init();

      //!identity
      static QSDArray<2> id;

      //!parity/sign
      static QSDArray<2> s;

      //!creator of an up-spin particle
      static QSDArray<2> cu;

      //!creator of a down-spin particle
      static QSDArray<2> cd;

      //!creator of an up-spin particle: with sign
      static QSDArray<2> cus;

      //!creator of a down-spin particle: with sign
      static QSDArray<2> cds;

      //!annihilator of an up-spin particle
      static QSDArray<2> au;

      //!annihilator of a down-spin particle
      static QSDArray<2> ad;

      //!annihilator of an up-spin particle: with sign
      static QSDArray<2> aus;

      //!annihilator of a down-spin particle: with sign
      static QSDArray<2> ads;

      //!creator of an up-down pair
      static QSDArray<2> cucd;

      //!creator of an up annihilator of an up
      static QSDArray<2> cuau;

      //!creator of an up annihilator of a down
      static QSDArray<2> cuad;

      //!creator of a down annihilator of an up
      static QSDArray<2> cdau;

      //!creator of a down annihilator of an down
      static QSDArray<2> cdad;

      //!annihilator of an up annihilator of a down
      static QSDArray<2> auad;

      //!triple operator with signature of an up spin creator
      static QSDArray<2> tcu;

      //!triple operator with signature of a down spin creator
      static QSDArray<2> tcd;

      //!triple operator with signature of an down spin annihilator
      static QSDArray<2> tad;

      //!triple operator with signature of an up spin annihilator
      static QSDArray<2> tau;

      //!local two-body operator
      static QSDArray<2> ltb;

      //!local one-body operator
      static QSDArray<2> lob;

      //!print out all the operators
      static void print();

};

#endif
