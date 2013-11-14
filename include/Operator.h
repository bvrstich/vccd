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

      static void init(const DArray<2> &,const DArray<4> &);

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

      //!complementary triple operator with signature of an up spin creator
      static std::vector< std::vector< QSDArray<2> > > tcuf;

      //!complementary triple operator with signature of a down spin creator
      static std::vector< std::vector< QSDArray<2> > > tcdf;

      //!complementary triple operator with signature of an down spin annihilator
      static std::vector< std::vector< QSDArray<2> > > tadf;

      //!complementary triple operator with signature of an up spin annihilator
      static std::vector< std::vector< QSDArray<2> > > tauf;

      //!local operator
      static std::vector< QSDArray<2> > local;

      //!print out all the operators
      static void print();

};

#endif
