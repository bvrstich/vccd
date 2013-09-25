#ifndef T2_2_MPO_H
#define T2_2_MPO_H

#include <iostream>
#include <cstdlib>

using std::ostream;

class Ostate;

/**
 * @author Brecht Verstichel
 * @date 24-08-2013
 * This class just collects the information needed for incoming and outgoing states of a MPO.
 */
class T2_2_mpo {

   friend ostream &operator<<(ostream &output,const T2_2_mpo &list_p);

   public:

      //constructor
      T2_2_mpo(int,int);

      //copy constructor
      T2_2_mpo(const T2_2_mpo &);

      //destructor
      virtual ~T2_2_mpo();

      void push_anni_down_anni_up(int,int,int,int,int);

      void push_crea_up_crea_down(int,int,int,int,int);

      void push_anni_up(int,int,int,int,int);

      void push_anni_down_s(int,int,int,int,int);

      void push_crea_up_s(int,int,int,int,int);

      void push_crea_down(int,int,int,int,int);

      void push_id(int,int,int,int,int);

      void push_single_in_complement(int,const Ostate &,int,const Ostate &,int);

      void push_single_out_complement(int,const Ostate &,int,const Ostate &,int);

      void push_double_complement(int,const Ostate &,int,const Ostate &,int);

      template<class Q>
         double get(const MPO<Q> &,int,int,int,int) const;

   private:

      std::vector< std::vector<int> > *list;

      int **ij2o;
      std::vector< std::vector<int> > o2ij;

      int **ab2v;
      std::vector< std::vector<int> > v2ab;

      int **ov2s;
      std::vector< std::vector<int> > s2ov;

      int no;
      int nv;

};

#endif
