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
class eMPS : public std::vector< MPS<Quantum> > {

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

      void mult_exc(const MPO<Quantum> &O,const MPS<Quantum> &hf);

      eMPS &mult_gen(const MPO<Quantum> &O,const MPS<Quantum> &hf,const eMPS &ccd);

      double inprod(const MPO<Quantum> &H,const MPS<Quantum> &hf,const eMPS &) const;

      double dot(const eMPS &) const;

      MPS<Quantum> expand(const MPS<Quantum> &,int,int D) const;

      void print_dim(int);

      void print_tot_dim(int);

   private:

      //!the order of all the elements
      std::vector<int> order;

      //!the order of the lowest order term
      int ord;

      //!list containing the dimensions contributing to the different orders in the expansion
      std::vector<int> cutoff;

};

#endif
