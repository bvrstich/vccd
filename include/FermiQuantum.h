#ifndef SPINQUANTUM_H
#define SPINQUANTUM_H

#include <iostream>
#include <iomanip>

/**
 * class which defines the quantumnumbers of the system, this is the case of spinless fermions
 */
class FermiQuantum
{

   public:

      /**
       * empty constructor, set particle number to zero
       */
      FermiQuantum(){

         N = 0;

      }

      /**
       * constructor with input, set particle number
       * @param N_i input quantumnumber
       */
      FermiQuantum(int N_i){

         N = N_i;

      }

      /**
       * copy constructor
       */
      FermiQuantum(const FermiQuantum &qn_copy){

         N = qn_copy.gN();

      }

      /**
       * @return the quantumnumber
       */
      int gN() const {

         return N;

      }

      /**
       *  equality operator overloaded
       * @param qn_i input
       * @return true if input == *this
       */
      inline bool operator==(const FermiQuantum &qn_i) const { 

         return (N == qn_i.gN());

      }

      /**
       * inequality operator overloaded
       * @param qn_i input
       * @return true if input != *this
       */
      inline bool operator!=(const FermiQuantum& qn_i) const {

         return (N != qn_i.gN());

      }

      /**
       * < comparison operator overloaded
       * @param qn_i input
       * @return true if *this < input
       */
      inline bool operator<(const FermiQuantum& qn_i) const {

         return (N < qn_i.gN());

      }

      /**
       * > comparison operator overloaded
       * @param qn_i input
       * @return true if *this > input
       */
      inline bool operator>(const FermiQuantum& qn_i) const { 

         return (N > qn_i.gN());

      }

      /**
       * operator acting on quantumnumbers
       * @param qn_i input
       * @return new FermiQuantum object with N = *this + input
       */
      inline FermiQuantum operator*(const FermiQuantum& qn_i) const {

         return FermiQuantum(N + qn_i.gN());

      }

      /**
       * overload the + operator: basically makes a copy of the input object
       * @param input object q
       */
      friend FermiQuantum operator+ (const FermiQuantum& q) {

         return FermiQuantum(q.gN()); 

      }

      /**
       * overload the - operator: basically makes a copy of the input object with negative sign
       * @param input object q
       */
      friend FermiQuantum operator-(const FermiQuantum& q) { 

         return FermiQuantum(-q.gN()); 

      }

      /**
       * overload output stream operator
       */
      friend std::ostream& operator<< (std::ostream& ost, const FermiQuantum& q) {

         ost << "(" << std::setw(2) << q.N << ")";

         return ost;

      }

      inline bool parity() const { 

         return N & 1; 
         
      }


      /**
       * @return a FermiQuantum object initialized on zero
       */
      const static FermiQuantum zero() {

         return FermiQuantum(0); 

      }

   private:

      //! the number of particles
      int N;

};

#endif // SPINQUANTUM_H
