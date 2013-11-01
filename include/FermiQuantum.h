#ifndef FERMIQUANTUM_H
#define FERMIQUANTUM_H

#include <iostream>
#include <iomanip>
#include <boost/serialization/serialization.hpp>

/**
 * class which defines the quantumnumbers of the system, this is the case of spin-1/2 fermions
 */
class FermiQuantum
{

   public:

      /**
       * empty constructor, set particle number to zero
       */
      FermiQuantum(){

         n_up = 0;
         n_down = 0;

      }

      /**
       * constructor with input, set particle number
       * @param n_i input quantumnumber
       */
      FermiQuantum(int n_up_i,int n_down_i){

         n_up = n_up_i;
         n_down = n_down_i;

      }

      /**
       * copy constructor
       */
      FermiQuantum(const FermiQuantum &qn_copy){

         n_up = qn_copy.gn_up();
         n_down = qn_copy.gn_down();

      }

      /**
       * @return the number of spin up particles
       */
      int gn_up() const {

         return n_up;

      }

      /**
       * @return the number of spin down particles
       */
      int gn_down() const {

         return n_down;

      }

      /**
       *  equality operator overloaded
       * @param qn_i input
       * @return true if input == *this
       */
      inline bool operator==(const FermiQuantum &qn_i) const { 

         return ( (n_up == qn_i.gn_up()) && (n_down == qn_i.gn_down()) );

      }

      /**
       * inequality operator overloaded
       * @param qn_i input
       * @return true if input != *this
       */
      inline bool operator!=(const FermiQuantum& qn_i) const {

         return ( (n_up != qn_i.gn_up()) || (n_down != qn_i.gn_down()) );

      }

      /**
       * < comparison operator overloaded
       * @param qn_i input
       * @return true if *this < input
       */
      inline bool operator<(const FermiQuantum& qn_i) const {

         if(n_up < qn_i.gn_up())
            return true;
         else{

            if(n_up == qn_i.gn_up()){

               if(n_down < qn_i.gn_down())
                  return true;
               else
                  return false;

            }
            else
               return false;

         }

      }

      /**
       * > comparison operator overloaded
       * @param qn_i input
       * @return true if *this > input
       */
      inline bool operator>(const FermiQuantum& qn_i) const { 

         if(n_up > qn_i.gn_up())
            return true;
         else{

            if(n_up == qn_i.gn_up()){

               if(n_down > qn_i.gn_down())
                  return true;
               else
                  return false;

            }
            else
               return false;

         }

      }

      /**
       * operator acting on quantumnumbers
       * @param qn_i input
       * @return new FermiQuantum object with n = *this + input
       */
      inline FermiQuantum operator*(const FermiQuantum& qn_i) const {

         return FermiQuantum(n_up + qn_i.gn_up(),n_down + qn_i.gn_down());

      }

      /**
       * overload the + operator: basically makes a copy of the input object
       * @param input object q
       */
      friend FermiQuantum operator+ (const FermiQuantum& q) {

         return FermiQuantum(q.gn_up(),q.gn_down()); 

      }

      /**
       * overload the - operator: basically makes a copy of the input object with negative sign
       * @param input object q
       */
      friend FermiQuantum operator-(const FermiQuantum& q) { 

         return FermiQuantum(-q.gn_up(),-q.gn_down()); 

      }

      /**
       * unused function that has to be defined
       */
      inline bool parity() const { 

         return true; 
         
      }


      /**
       * output stream operator overloaded
       */
      friend std::ostream& operator<< (std::ostream& ost, const FermiQuantum& q) {

         ost << "(" << std::setw(2) << q.gn_up() << "," << std::setw(2) << q.gn_down() << ")";

         return ost;

      }

      /**
       * @return a FermiQuantum object initialized on zero
       */
      const static FermiQuantum zero() {

         return FermiQuantum(0,0); 

      }

      /**
       *  create an up particle to the quantumnumber (actually means -1 for outgoing)
       */
      void crea_up(){
         
         n_up -= 1;

      }

      /**
       * crea_down an up particle to the quantumnumber (actually means -1 for outgoing)
       */
      void crea_down(){
         
         n_down -= 1;

      }

      /**
       * anni_up an up particle to the quantumnumber (actually means +1 for outgoing)
       */
      void anni_up(){
         
         n_up += 1;

      }

      /**
       * anni_down an down particle to the quantumnumber (actually means +1 for outgoing)
       */
      void anni_down(){
         
         n_down += 1;

      }

   private:

      friend class boost::serialization::access;

      template <class Archive>
         void serialize(Archive& ar, const unsigned int version) {
            
            ar & n_up & n_down;
         
         }

      //! the number of up-particles
      int n_up;

      //! the number of down-particles
      int n_down;

};

#endif // 
