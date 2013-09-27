#ifndef EMPS_H
#define EMPS_H

#include <iostream>
#include <cstdlib>

using std::ostream;

/**
 * @author Brecht Verstichel
 * @date 26-09-2013
 * This is a class written for making dealing with exponential MPO/MPS stuff more efficient 
 */
class eMPS {

   friend ostream &operator<<(ostream &output,const eMPS &emps_p);

   public:

      //constructor
      eMPS();

      eMPS(const MPO<Quantum> &,const MPS<Quantum> &,const std::vector<int> &);

      //copy constructor
      eMPS(const eMPS &);

      //destructor
      virtual ~eMPS();

      void update(const MPO<Quantum> &T,const MPS<Quantum> &hf);

      const std::vector<int> &gcutoff() const;

      double eval(const MPO<Quantum> &H,const MPS<Quantum> &hf) const;

      double norm() const;

      const std::vector<int> &gorder() const;

      void mult(const MPO<Quantum> &O,const MPS<Quantum> &hf);

      double inprod(const MPO<Quantum> &H,const MPS<Quantum> &hf,const eMPS &) const;

      double dot(const eMPS &) const;

      MPS<Quantum> expand(const MPS<Quantum> &,int D = 0) const;

   private:

      //!vector containing the higher order terms: T*hf, T^2*hf, etc...
      std::vector< MPS<Quantum> > term;

      //!the order of all the elements
      std::vector<int> order;

      //!the order of the lowest order term
      int ord;

      //!list containing the dimensions contributing to the different orders in the expansion
      std::vector<int> cutoff;

};

#endif
