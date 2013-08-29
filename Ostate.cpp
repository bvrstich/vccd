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

/**
 * get the complementary operator between in and out, when a pair is coming and one is going out and give the correct sign to rearrange to normal order
 */
std::vector<int> Ostate::get_single_complement(int site,const Ostate &in,const Ostate &out,const DArray<4> &V,double &Veff){

   //lets call in i,j and out l
   int i = oplist[in[1]][0];
   int j = oplist[in[0]][0];

   int si = oplist[in[1]][1];
   int sj = oplist[in[0]][1];

   int ai = oplist[in[1]][2];
   int aj = oplist[in[0]][2];

   int l = oplist[out[0]][0];
   int sl = oplist[out[0]][1];
   int al = oplist[out[0]][2];

   vector<int> comp;

   if(si == 0 && ai == 0){//first site create up spin

      if(sj == 0 && aj == 0){//second site create up spin

         if(sl == 0 || al == 1){//only coupling if annihilator with spin  up

            comp.resize(2);
            comp[0] = 0;
            comp[1] = 1;

            Veff = V(i,j,site,l) - V(i,j,l,site);

            return comp;

         }
         else{

            Veff = 0.0;
            return comp;

         }

      }
      else if(sj == 0 && aj == 1){//second site annihilate up spin

         if(sl == 0 && al == 0){

            //put in anni up
            comp.resize(2);
            comp[0] = 0;
            comp[1] = 1;

            Veff = V(i,l,j,site) - V(i,l,site,j);

            return comp;

         }
         else if(sl == 0 && al == 1){

            //put in crea up
            comp.resize(2);
            comp[0] = 0;
            comp[1] = 0;

            Veff = V(i,site,j,l) - V(i,site,l,j);

            return comp;

         }
         else if(sl == 1 && al == 0){//create down

            //put in anni down
            comp.resize(2);
            comp[0] = 1;
            comp[1] = 1;

            Veff = V(i,l,j,site);

            return comp;

         }
         else if(sl == 1 && al == 1){//annihilate down

            //put in crea down
            comp.resize(2);
            comp[0] = 1;
            comp[1] = 0;

            Veff = V(i,site,j,l);

            return comp;

         }

      }
      else if(sj == 1 && aj == 0){//second site create down spin

         if(al == 0){//impossible

            Veff = 0.0;
            return comp;

         }
         else{

            if(sl == 0){//anni up

               //put in anni down
               comp.resize(2);
               comp[0] = 1;
               comp[1] = 1;

               Veff = -V(i,j,l,site);

               return comp;

            }
            else{//anni down

               //put in anni up
               comp.resize(2);
               comp[0] = 0;
               comp[1] = 1;

               Veff = V(i,j,site,l);

               return comp;

            }

         }

      }
      else{//second site anni down spin

         if(al == 0 && sl == 1){//out is crea down

            //put in anni up
            comp.resize(2);
            comp[0] = 0;
            comp[1] = 1;

            Veff = -V(i,l,site,j);

            return comp;


         }
         else if(al == 1 && sl == 0){

            //put in anni up
            comp.resize(2);
            comp[0] = 0;
            comp[1] = 1;

            Veff = -V(i,site,l,j);

            return comp;

         }
         else{

            Veff = 0.0;
            return comp;

         }

      }

   }//end first site
   else if(si == 1 && ai == 0){//first site create down spin

      if(sj == 0 && aj == 0){//second site create up spin

         if(al == 1){//only coupling if annihilator 

            if(sl == 0){//if spin up

               comp.resize(2);
               comp[0] = 0;
               comp[1] = 1;

               Veff = V(i,j,site,l);

               return comp;

            }
            else{//if spin down output

               comp.resize(2);
               comp[0] = 0;
               comp[1] = 1;

               Veff = -V(i,j,l,site);

               return comp;

            }

         }
         else{

            Veff = 0.0;
            return comp;

         }

      }
      else if(sj == 0 && aj == 1){//second site annihilate up spin

         if(sl == 0 && al == 0){//out create up spin

            comp.resize(2);
            comp[0] = 1;
            comp[1] = 1;

            Veff = V(l,i,j,site);

            return comp;

         }
         else if(sl == 1 && al == 1){

            comp.resize(2);
            comp[0] = 0;
            comp[1] = 0;

            Veff = -V(site,i,j,l);

            return comp;

         }
         else{

            Veff = 0.0;
            return comp;

         }

      }
      else if(sj == 1 && aj == 0){//second site create down spin

         if(al == 1 && sl == 1){//only possibility

            //put in anni down
            comp.resize(2);
            comp[0] = 1;
            comp[1] = 1;

            Veff = V(i,j,site,l) - V(i,j,l,site);

            return comp;

         }
         else{

            Veff = 0.0;
            return comp;

         }

      }
      else{//second site anni down spin

         if(al == 1 && sl == 1){

            //put in crea down
            comp.resize(2);
            comp[0] = 1;
            comp[1] = 0;

            Veff = V(i,site,j,l) - V(i,site,l,j);

            return comp;

         }
         else if(al == 0 && sl == 1){

            //put in anni down
            comp.resize(2);
            comp[0] = 1;
            comp[1] = 1;

            Veff = V(i,l,j,site) - V(i,l,site,j);

            return comp;

         }
         else if(al == 1 && sl == 0){//anni up

            //put in crea up 
            comp.resize(2);
            comp[0] = 0;
            comp[1] = 0;

            Veff = V(site,i,l,j);

            return comp;

         }
         else if(al == 1 && sl == 0){//crea up

            //put in anni up 
            comp.resize(2);
            comp[0] = 0;
            comp[1] = 1;

            Veff = V(l,i,site,j);

            return comp;

         }

      }

   }//end first site
   else if(si == 0 && ai == 1){//first site anni up spin

      if(sj == 0 && aj == 0){//second site create up spin

         if(al == 1 && sl == 0){//out anni up spin

            comp.resize(2);
            comp[0] = 0;
            comp[1] = 0;

            Veff = V(j,site,i,l) - V(j,site,l,i);

            return comp;

         }
         else if(al == 0 && sl == 0){//out crea up spin

            comp.resize(2);
            comp[0] = 0;
            comp[1] = 1;

            Veff = V(j,l,i,site) - V(j,l,site,i);

            return comp;

         }
         else if(al == 1 && sl == 1){//anni down

            comp.resize(2);
            comp[0] = 1;
            comp[1] = 0;

            Veff = V(j,site,i,l);

            return comp;

         }
         else if(al == 0 && sl == 1){//crea down

            comp.resize(2);
            comp[0] = 1;
            comp[1] = 1;

            Veff = V(j,l,i,site);

            return comp;

         }

      }
      else if(sj == 0 && aj == 1){//second site annihilate up spin

         if(sl == 0 && al == 0){//only possib

            comp.resize(2);
            comp[0] = 0;
            comp[1] = 0;

            Veff = V(site,l,i,j) - V(site,l,j,i);

            return comp;

         }
         else{

            Veff = 0.0;
            return comp;

         }

      }
      else if(sj == 1 && aj == 0){//second site create down spin

         if(al == 0 && sl == 0){

            comp.resize(2);
            comp[0] = 1;
            comp[1] = 1;

            Veff = -V(l,j,i,site);

            return comp;

         }
         else if(al == 1 && sl == 1){

            comp.resize(2);
            comp[0] = 0;
            comp[1] = 0;

            Veff = -V(site,j,i,l);

            return comp;

         }
         else{

            Veff = 0.0;
            return comp;

         }

      }
      else{//second site anni down spin

         if(al == 0){

            if(sl == 0){

               comp.resize(2);
               comp[0] = 1;
               comp[1] = 0;

               Veff = -V(i,j,l,site);

               return comp;

            }
            else{

               comp.resize(2);
               comp[0] = 0;
               comp[1] = 0;

               Veff = V(i,j,site,l);

               return comp;

            }

         }
         else{

            Veff = 0.0;
            return comp;

         }

      }

   }
   else{//first site anni down spin

      if(sj == 0 && aj == 0){//second site create up spin

         if(al == 0 && sl == 1){//out crea down spin

            comp.resize(2);
            comp[0] = 0;
            comp[1] = 1;

            Veff = -V(j,l,site,i);

            return comp;

         }
         else if(al == 1 && sl == 0){//out anni up spin

            comp.resize(2);
            comp[0] = 1;
            comp[1] = 0;

            Veff = -V(j,site,l,i);

            return comp;

         }
         else{

            Veff = 0.0;
            return comp;

         }

      }
      else if(sj == 0 && aj == 1){//second site annihilate up spin

         if(al == 0){//only creators

            if(sl == 0){

               comp.resize(2);
               comp[0] = 1;
               comp[1] = 0;

               Veff = V(i,j,site,l);

               return comp;

            }
            else{

               comp.resize(2);
               comp[0] = 0;
               comp[1] = 0;

               Veff = -V(i,j,l,site);

               return comp;

            }

         }
         else{

            Veff = 0.0;
            return comp;

         }

      }
      else if(sj == 1 && aj == 0){//second site create down spin

         if(al == 0 && sl == 1){

            comp.resize(2);
            comp[0] = 1;
            comp[1] = 1;

            Veff = V(j,l,i,site) - V(j,l,site,i);

            return comp;

         }
         else if(al == 1 && sl == 1){

            comp.resize(2);
            comp[0] = 1;
            comp[1] = 0;

            Veff = V(j,site,i,l) - V(j,site,l,i);

            return comp;

         }
         else if(al == 0 && sl == 0){//out is create up

            comp.resize(2);
            comp[0] = 0;
            comp[1] = 1;

            Veff = V(i,site,j,l);

            return comp;

         }
         else{//out is anni up

            comp.resize(2);
            comp[0] = 0;
            comp[1] = 0;

            Veff = V(site,j,i,l);

            return comp;

         }

      }
      else{//second site anni down spin

         if(al == 0 && sl == 1){

            comp.resize(2);
            comp[0] = 1;
            comp[1] = 0;

            Veff = V(site,l,i,j) - V(site,l,j,i);

            return comp;

         }
         else{

            Veff = 0.0;
            return comp;

         }

      }

   }

}

/**
 * get the complementary pair operator between in and out, when one operator is coming and one is going out.
 * @return a vector containing the info about which operator is the complement:
 * if vector size is zero, no complement
 * if vector size is one: v[0] = 0 -> a_down a_up
 *                        v[0] = 1 -> a^+_up a^+_down
 *                        v[0] = 2 -> a^+_up a_down
 *                        v[0] = 3 -> a^+_down a_up
 * if vector size is 2: operator is sum of up up and down down:
 * coefficients and signs are returned in Veff
 */
std::vector<int> Ostate::get_double_complement(int site,const Ostate &in,const Ostate &out,const DArray<4> &V,std::vector<double> &Veff){

   //lets call in i and out l
   int i = oplist[in[0]][0];
   int si = oplist[in[0]][1];
   int ai = oplist[in[0]][2];

   int l = oplist[out[0]][0];
   int sl = oplist[out[0]][1];
   int al = oplist[out[0]][2];

   vector<int> comp;

   if(si == 0 && ai == 0){//in create up spin

      if(sl == 0 && al == 0)//out create up spin: not possible
         return comp;
      else if(sl == 1 && al == 0){//out create down spin

         comp.push_back(0);
         Veff[0] = V(i,l,site,site);

         return comp;

      }
      else if(sl == 0 && al == 1){//out annihilate up spin

         comp.resize(2);
         Veff[0] = V(i,site,site,l) - V(i,site,l,site);
         Veff[1] = -V(i,site,l,site);

         return comp;

      }
      else{//out anni down spin

         comp.push_back(3);
         Veff[0] = V(i,site,site,l);

         return comp;

      }

   }
   else if(si == 1 && ai == 0){//in create down spin

      if(sl == 0 && al == 0){//out create up spin

         comp.push_back(0);
         Veff[0] = -V(l,i,site,site);

         return comp;

      }
      else if(sl == 1 && al == 0)//out create down spin: no possible
         return comp;
      else if(sl == 0 && al == 1){//out annihilate up spin

         comp.push_back(2);
         Veff[0] = V(site,i,l,site);

         return comp;

      }
      else{//out anni down spin

         comp.resize(2);
         Veff[0] = -V(site,i,site,l);
         Veff[1] = V(i,site,site,l) - V(i,site,l,site);

         return comp;

      }

   }
   else if(si == 0 && ai == 1){//anni up spin

      if(sl == 0 && al == 0){//out create up spin

         comp.resize(2);
         Veff[0] = V(i,site,site,l) - V(i,site,l,site);
         Veff[1] = -V(l,site,i,site);

         return comp;

      }
      else if(sl == 1 && al == 0){//out create down spin:

         comp.push_back(2);
         Veff[0] = V(site,l,i,site);

         return comp;

      }
      else if(sl == 0 && al == 1)//out annihilate up spin: no possible
         return comp;
      else{//out anni down spin

         comp.push_back(1);
         Veff[0] = V(site,site,i,l);

         return comp;

      }

   }
   else{//in anni down spin

      if(sl == 0 && al == 0){//out create up spin

         comp.push_back(3);
         Veff[0] = V(l,site,site,i);

         return comp;

      }
      else if(sl == 1 && al == 0){//out create down spin:

         comp.resize(2);
         Veff[0] = -V(l,site,i,site);
         Veff[1] = V(i,site,site,l) - V(i,site,l,site);

         return comp;

      }
      else if(sl == 0 && al == 1){//out annihilate up spin: no possible

         comp.push_back(1);
         Veff[0] = V(site,site,i,l);

         return comp;

      }
      else//out anni down spin: no possible
         return comp;

   }

}

/**
 * for an incoming pair of operators in 'in', get the complement pair that closes the expression down
 * if vector size is zero, no complement
 * if vector size is one: v[0] = 0 -> a_down a_up
 *                        v[0] = 1 -> a^+_up a^+_down
 *                        v[0] = 2 -> a^+_up a_down
 *                        v[0] = 3 -> a^+_down a_up
 * if vector size is 2: operator is sum of up up and down down:
 */
std::vector<int> Ostate::get_closing_pair(int site,const Ostate &in,const DArray<4> &V,vector<double> &val){

   //lets call in i,j
   int i = oplist[in[1]][0];
   int j = oplist[in[0]][0];

   int si = oplist[in[1]][1];
   int sj = oplist[in[0]][1];

   int ai = oplist[in[1]][2];
   int aj = oplist[in[0]][2];

   vector<int> comp;

   if(si == 0 && ai == 0){//create up

      if(sj == 0 && aj == 0)//create up
         return comp;
      else if(sj == 1 && aj == 0){//create down

         comp.push_back(0);
         val[0] = V(i,j,site,site);

         return comp;

      }
      else if(sj == 0 && aj == 1){//annihilate up

         comp.resize(2);
         val[0] = V(i,site,j,site) - V(i,site,site,j);
         val[1] = V(i,site,j,site);

         return comp;

      }
      else{

         comp.push_back(3);
         val[0] = -V(i,site,site,j);

         return comp;

      }

   }
   else if(si == 1 && ai == 0){//create down

      if(sj == 0 && aj == 0){//create up

         comp.push_back(0);
         val[0] = -V(j,i,site,site);

         return comp;

      }
      else if(sj == 1 && aj == 0)//create down
         return comp;
      else if(sj == 0 && aj == 1){//annihilate up

         comp.push_back(2);
         val[0] = -V(site,i,site,j);

         return comp;

      }
      else{

         comp.resize(2);
         val[0] = V(site,i,site,j);
         val[1] = V(i,site,j,site) - V(i,site,site,j);

         return comp;

      }

   }
   else if(si == 0 && ai == 1){//annihilate up

      if(sj == 0 && aj == 0){//create up

         comp.resize(2);
         val[0] = V(j,site,i,site) - V(j,site,site,i);
         val[1] = V(j,site,i,site);

         return comp;

      }
      else if(sj == 1 && aj == 0){//create down

         comp.push_back(2);
         val[0] = V(site,j,i,site);

         return comp;

      }
      else if(sj == 0 && aj == 1)//annihilate up
         return comp;
      else{

         comp.push_back(1);
         val[0] = V(j,i,site,site);

         return comp;

      }

   }
   else{//annihilate down

      if(sj == 0 && aj == 0){//create up

         comp.push_back(3);
         val[0] = V(j,site,i,site);

         return comp;

      }
      else if(sj == 1 && aj == 0){//create down

         comp.resize(2);
         val[0] = V(site,j,site,i);
         val[1] = V(j,site,i,site) - V(j,site,site,i);

         return comp;

      }
      else if(sj == 0 && aj == 1){//annihilate up

         comp.push_back(1);
         val[0] = -V(site,site,i,j);

         return comp;

      }
      else 
         return comp;

   }

}

/**
 * for a one body operator, get the correct closing operator and return the value and sign
 */
std::vector<int> Ostate::get_closing_single(int site,const Ostate &in,const DArray<2> &t,double &val){

   //lets call in i,j
   int i = oplist[in[0]][0];
   int si = oplist[in[0]][1];
   int ai = oplist[in[0]][2];

   vector<int> comp;

   if(si == 0 && ai == 0){//create up

      comp.resize(2);
      comp[0] = 0;//up
      comp[1] = 1;//anni

      val = t(i,site);

      return comp;

   }
   else if(si == 1 && ai == 0){//create down

      comp.resize(2);
      comp[0] = 1;//down
      comp[1] = 1;//anni

      val = t(i,site);

      return comp;

   }
   else if(si == 0 && ai == 1){//anni up

      comp.resize(2);
      comp[0] = 0;//up
      comp[1] = 0;//crea

      val = t(site,i);

      return comp;

   }
   else{//anni down

      comp.resize(2);
      comp[0] = 1;//down
      comp[1] = 0;//crea

      val = t(site,i);

      return comp;

   }

}
