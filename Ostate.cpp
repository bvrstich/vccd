#include <iostream>
#include <fstream>
#include <cmath>

using std::endl;
using std::ostream;
using std::ofstream;
using std::ifstream;

#include "include.h"

std::vector< std::vector<int> > Ostate::oplist;

/**
 * static fucntion which constructs the list relating the operator index to the actual operator:
 * 0 == id
 * 1 -> L = create up on site
 * L+1 -> 2L = create down on site 
 * 2L + 1 -> 3L = anni up on site 
 * 3L + 1 -> 4L = anni down on site
 * @param L length of the chain 
 */
void Ostate::construct_oplist(int L){

   std::vector<int> v(1);//Id
   v[0] = 0;

   oplist.push_back(v);

   v.resize(3);//Id

   for(int i = 0;i < L;++i){

      v[0] = i;//site
      v[1] = 0;//spin up
      v[2] = 0;//create

      oplist.push_back(v);

   }

   for(int i = 0;i < L;++i){

      v[0] = i;//site
      v[1] = 1;//spin down
      v[2] = 0;//create

      oplist.push_back(v);

   }

   for(int i = 0;i < L;++i){

      v[0] = i;//site
      v[1] = 0;//spin up
      v[2] = 1;//anni

      oplist.push_back(v);

   }

   for(int i = 0;i < L;++i){

      v[0] = i;//site
      v[1] = 1;//spin down
      v[2] = 1;//anni

      oplist.push_back(v);

   }

}

/**
 * static function: prints the oplist
 */
void Ostate::print_oplist(){

   for(int i = 0;i < oplist.size();++i)
      cout << i << "\t" << oplist[i] << endl;

}

/**
 * constructor
 * @param number of sites
 */
Ostate::Ostate(){ }

/**
 * constructor 
 * @param L number of sites
 * @param nop number of operators on the site
 */
Ostate::Ostate(int nop){

   op.resize(nop);

}

/**
 * copy constructor 
 * @param ostate_c object to be copied
 */
Ostate::Ostate(const Ostate &ostate_c){

   op = ostate_c.op;

}

Ostate::~Ostate(){ }

/**
 * overload equality operator
 */
Ostate &Ostate::operator=(const Ostate &ostate_c){

   op = ostate_c.op;

}

ostream &operator<<(ostream &output,const Ostate &ostate_p){

   for(int i = 0;i < ostate_p.size();++i)
      output << oplist[ostate_p.op[i]] << endl;

}

/**
 * add the identity operator
 */
void Ostate::push_id(){

   op.

}
