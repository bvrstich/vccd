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
class Ostate{

   friend ostream &operator<<(ostream &output,const Ostate &ostate_p);

   public:

      //constructor
      Ostate(int,bool,bool);

      //copy constructor
      Ostate(const Ostate &);

      //destructor
      virtual ~Ostate();

      //overload equality operator
      Ostate &operator=(const Ostate &);

      int gsite() const;

      bool gspin() const;

      bool gact() const;

   private:
      
      //!the site on which the operator acts
      int site;

      //!up (true) or down (false) spin
      bool spin;

      //!creation (true) or annihilation (false)
      bool act;

};

#endif
