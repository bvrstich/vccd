#ifndef E_EMPS_H
#define E_EMPS_H

#include <iostream>
#include <cstdlib>

using std::ostream;

class eMPS;

/**
 * @author Brecht Verstichel
 * @date 26-09-2013
 * This is a class written for making dealing with exponential MPO/MPS stuff more efficient. What happens if you multiply an exponential of an exponential ...
 */
class e_eMPS : public std::vector< eMPS > {

   friend ostream &operator<<(ostream &output,const e_eMPS &lsemps_p);

   public:

      e_eMPS(const MPO<Quantum> &O,const eMPS &,const MPS<Quantum> &);

      //copy constructor
      e_eMPS(const e_eMPS &);

      //destructor
      virtual ~e_eMPS();

      void fillE(DArray<2> &,const MPO<Quantum> &H,const MPS<Quantum> &hf,const eMPS &ccd) const;

      void fillN(DArray<2> &,const eMPS &ccd) const;

   private:

};

#endif
