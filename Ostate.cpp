#include <iostream>
#include <fstream>
#include <cmath>

using std::endl;
using std::ostream;
using std::ofstream;
using std::ifstream;

#include "include.h"

std::vector< std::vector<int> > Ostate::oplist;
int Ostate::L;

/**
 * static fucntion which constructs the list relating the operator index to the actual operator:
 * 0 == id
 * 1 -> L = create up on site
 * L+1 -> 2L = create down on site 
 * 2L + 1 -> 3L = anni up on site 
 * 3L + 1 -> 4L = anni down on site
 * @param L_in length of the chain 
 */
void Ostate::construct_oplist(int L_in){

   L = L_in;

   std::vector<int> v(1);//Id
   v[0] = 0;

   oplist.push_back(v);

   v.resize(3);

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
 */
Ostate::Ostate() : std::vector<int>(){ }

/**
 * constructor 
 * @param L number of sites
 * @param nop number of operators on the site
 */
Ostate::Ostate(int nop) : std::vector<int>(nop){ }

/**
 * copy constructor 
 * @param ostate_c object to be copied
 */
Ostate::Ostate(const Ostate &ostate_c) : std::vector<int>(ostate_c){ }

Ostate::~Ostate(){ }

ostream &operator<<(ostream &output,const Ostate &ostate_p){

   for(int i = 0;i < ostate_p.size();++i)
      output << Ostate::oplist[ostate_p[i]] << endl;

   return output;

}

/**
 * add the identity operator
 */
void Ostate::push_id(){

   push_back(0);

}

/**
 * add a creation operator of a spin up particle on site
 */
void Ostate::push_crea_up(int site){

   push_back(1 + site);

}

/**
 * add a creation operator of a spin down particle on site
 */
void Ostate::push_crea_down(int site){

   push_back(L + 1 + site);

}

/**
 * add a annihilation operator of a spin up particle on site
 */
void Ostate::push_anni_up(int site){

   push_back(2*L + 1 + site);

}

/**
 * add a annihilation operator of a spin up particle on site
 */
void Ostate::push_anni_down(int site){

   push_back(3*L + 1 + site);

}

/**
 * @return spin of the 'n'the operator: 0 is up, 1 is down
 */
int Ostate::gspin(int n) const{

   return oplist[(*this)[n]][1];

}

/**
 * @return site of the 'n'the operator
 */
int Ostate::gsite(int n) const{

   return oplist[(*this)[n]][0];

}

/**
 * @return action f the 'n'the operator: 0 is create, 1 is annihilate
 */
int Ostate::gact(int n) const{

   return oplist[(*this)[n]][2];

}

/**
 * @return true if the two operators are on the same site
 */
bool Ostate::is_pair() const {

   if(oplist[(*this)[0]][0] == oplist[(*this)[1]][0])
      return true;
   else
      return false;

}
