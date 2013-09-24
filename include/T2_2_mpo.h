#ifndef T2_2_MPO_H
#define T2_2_MPO_H

#include <iostream>
#include <cstdlib>

using std::ostream;

/**
 * @author Brecht Verstichel
 * @date 24-08-2013
 * This class just collects the information needed for incoming and outgoing states of a MPO.
 */
class T2_2_mpo : public std::vector<int> {

   friend ostream &operator<<(ostream &output,const T2_2_mpo &list_p);

   public:

      //constructor
      T2_2_mpo(int,int);

      //copy constructor
      T2_2_mpo(const T2_2_mpo &);

      //destructor
      virtual ~T2_2_mpo();

   private:

      std::vector< vector<int> > *list;

      int **ij2o;
      vector< vector<int> > o2ij;

      int **ab2v;
      vector< vector<int> > v2ab;

      int **ov2s;
      vector< vector<int> > s2ov;

      int no;
      int nv;

};

#endif
