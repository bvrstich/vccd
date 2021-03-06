#ifndef OSTATE_H
#define OSTATE_H

#include <iostream>
#include <cstdlib>

using std::ostream;

/**
 * @author Brecht Verstichel
 * @date 24-08-2013
 * This class just collects the information needed for incoming and outgoing states of a MPO.
 */
class Ostate : public std::vector<int> {

   friend ostream &operator<<(ostream &output,const Ostate &ostate_p);

   public:

      //constructor
      Ostate();

      Ostate(int);

      //copy constructor
      Ostate(const Ostate &);

      //destructor
      virtual ~Ostate();

      void resize(int);

      void push_id();

      void push_crea_up(int);

      void push_crea_down(int);

      void push_anni_up(int);

      void push_anni_down(int);

      int gspin(int) const;

      int gsite(int) const;

      int gact(int) const;

      bool is_pair() const;

      static std::vector<int> get_single_complement(int site,const Ostate &in,const Ostate &out,const DArray<4> &V,double &);

      static int get_single_complement_T2(int site,const Ostate &in,const Ostate &out,const DArray<4> &t,double &);

      static int get_single_complement_T2_bis(int site,const Ostate &in,const Ostate &out,const DArray<4> &t,double &);

      static int get_double_complement_T2(const Ostate &in,const Ostate &out,const DArray<4> &t,double &);

      static int transfer_pair_in_pair_out(const Ostate &in,const Ostate &out,const DArray<4> &V,double &);

      static std::vector<int> get_double_complement(int site,const Ostate &in,const Ostate &out,const DArray<4> &V,std::vector<double> &);
      
      static std::vector<int> get_closing_pair_in(int site,const Ostate &in,const DArray<4> &V,std::vector<double> &);

      static std::vector<int> get_closing_pair_out(int site,const Ostate &out,const DArray<4> &V,std::vector<double> &);

      static std::vector<int> get_closing_single(int site,const Ostate &in,const DArray<2> &t,double &);

      static void construct_oplist(int);

      static void print_oplist();

   private:

      //!vector containing the operators: 
      //!0 = id
      //!1 -> L = create up on site 
      //!L+1 -> 2L = create down on site 
      //!2L + 1 -> 3L = anni up on site 
      //!3L + 1 -> 4L = anni down
      static std::vector< std::vector<int> > oplist;

      //!number of sites
      static int L;
      
};

#endif
