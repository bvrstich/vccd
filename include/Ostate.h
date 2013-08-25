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

      using std::vector<int>::operator=;

      using std::vector<int>::operator[];

      void resize(int);

      void push_id();

      void push_crea_up(int);

      void push_crea_down(int);

      void push_anni_up(int);

      void push_anni_down(int);

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
