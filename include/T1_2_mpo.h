#ifndef T1_2_MPO_H
#define T1_2_MPO_H

#include <iostream>
#include <cstdlib>

using std::ostream;

/**
 * @author Brecht Verstichel
 * @date 24-08-2013
 * This class just collects the information needed for incoming and outgoing states of a MPO.
 */
class T1_2_mpo : public std::vector<int> {

   friend ostream &operator<<(ostream &output,const T1_2_mpo &ostate_p);

   public:

      //constructor
      T1_2_mpo(int,int);

      //copy constructor
      T1_2_mpo(const T1_2_mpo &);

      //destructor
      virtual ~T1_2_mpo();

   private:

      std::vector<int*> *list;

      int **ia2s;

      int **s2ia;

      int no;
      int nv;

};

#endif
