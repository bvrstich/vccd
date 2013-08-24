#include <iostream>
#include <fstream>
#include <cmath>

using std::endl;
using std::ostream;
using std::ofstream;
using std::ifstream;

#include "include.h"

/**
 * constructor 
 * @param site_i site of the state
 * @param spin_i spin of the state
 * @param act_i action on the state: create or annihilate
 */
Ostate::Ostate(int site_i,bool spin_i,bool act_i){

   this->site = site_i;
   this->spin = spin_i;
   this->act = act_i;

}

/**
 * copy constructor 
 * @param ostate_c object to be copied
 */
Ostate::Ostate(const Ostate &ostate_c){

   site = ostate_c.gsite();
   spin = ostate_c.gspin();
   act = ostate_c.gact();

}

/**
 * @return the site
 */
int Ostate::gsite(){

   return site;

}

/**
 * @return the spin
 */
bool Ostate::gspin(){

   return spin;

}

/**
 * @return the action on the site
 */
bool Ostate::gact(){

   return act;

}
