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

      static void clear();

   private:
      
      //!vector containing all the 4x4 operators on every site for incoming and outgoing.
      static std::vector< QSDArray<2> ** > op;

      //!vector containing the sparsity information
      static std::vector< bool* > sparse;

      //!number of terms on row and columns of MPO
      static std::vector< std::vector<int> > dim;

      //!number of sites
      static int L;

};

#endif
