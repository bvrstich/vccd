#include <iostream>
#include <fstream>
#include <cmath>

using std::endl;
using std::ostream;
using std::ofstream;
using std::ifstream;

#include "include.h"

/**
 * empty constructor
 */
eMPS::eMPS(){ }

/**
 * standard constructor, initializes the T and hf
 */
eMPS::eMPS(const MPO<Quantum> &T,const MPS<Quantum> &hf,const std::vector<int> &cutoff_in){

   cutoff = cutoff_in;

   term.resize(cutoff.size());

   //form the list of contributing terms in the expansion
   term[0] = T*hf;

   compress(term[0],mpsxx::Left,cutoff[0]);
   compress(term[0],mpsxx::Right,cutoff[0]);

   for(int i = 1;i < cutoff.size();++i){

      term[i] = T*term[i - 1];
      compress(term[i],mpsxx::Left,cutoff[i]);
      compress(term[i],mpsxx::Right,cutoff[i]);
      scal(1.0/(i + 1.0),term[i]);

   }

   order.resize(cutoff.size());

   for(int i = 0;i < order.size();++i)
      order[i] = i + 1;

   ord = 0;

}

/**
 * copy constructor
 */
eMPS::eMPS(const eMPS &emps){

   term = emps.term;

   cutoff = emps.gcutoff();

   order = emps.gorder();

   ord = emps.ord;

}

/**
 * destructor
 */
eMPS::~eMPS(){ }

/**
 * @return the list containing the contributions of the different dimensions to the expansion
 */
const std::vector<int> &eMPS::gcutoff() const {

   return cutoff;

}

/**
 * @return the list containing the contributions of the different dimensions to the expansion
 */
const std::vector<int> &eMPS::gorder() const {

   return order;

}


/**
 * evaluate the expectation value of an operator, assume the operator connects only 2 particle terms
 */
double eMPS::eval(const MPO<Quantum> &H,const MPS<Quantum> &hf) const {

   double E = mpsxx::inprod(mpsxx::Left,hf,H,hf) + 2.0 * mpsxx::inprod(mpsxx::Left,hf,H,term[0]);

   for(int i = 0;i < cutoff.size();++i){

      E += mpsxx::inprod(mpsxx::Left,term[i],H,term[i]);

      if(i + 1 < cutoff.size())
         E += 2.0 * mpsxx::inprod(mpsxx::Left,term[i],H,term[i + 1]);

   }

   return E;

}

/**
 * evaluate the expectation value of an operator, assume the operator connects only 2 particle terms
 */
double eMPS::norm() const {

   double nrm = 1.0;

   for(int i = 0;i < cutoff.size();++i)
       nrm += mpsxx::dot(mpsxx::Left,term[i],term[i]);

   return nrm;

}

/**
 * @return the exponential as one MPS
 */
MPS<Quantum> eMPS::expand(const MPS<Quantum> &hf,int D) const{

   MPS<Quantum> sum = term[0];

   //now sum all the terms together:
   for(int i = 1;i < cutoff.size();++i){

      axpy(1.0,term[i],sum);
      compress(sum,mpsxx::Left,D);
      compress(sum,mpsxx::Right,D);

   }

   axpy(1.0,hf,sum);

   return sum;

}

/**
 * update the eMPS for a new MPO form of the T operator
 */
void eMPS::update(const MPO<Quantum> &T,const MPS<Quantum> &hf){

   //form the list of contributing terms in the expansion
   term[0] = T*hf;

   compress(term[0],mpsxx::Left,cutoff[0]);
   compress(term[0],mpsxx::Right,cutoff[0]);

   for(int i = 1;i < cutoff.size();++i){

      term[i] = T*term[i - 1];
      compress(term[i],mpsxx::Left,cutoff[i]);
      compress(term[i],mpsxx::Right,cutoff[i]);
      scal(1.0/(i + 1.0),term[i]);

   }

}

/**
 * multiply this with an MPO
 * @param O input MPO
 */
void eMPS::mult(const MPO<Quantum> &O,const MPS<Quantum> &hf){

   if(ord == 0){

      for(int i = term.size() - 1;i > 0;--i){

         term[i] = O * term[i - 1];

         compress(term[i],mpsxx::Left,cutoff[i]);
         compress(term[i],mpsxx::Right,cutoff[i]);

      }

      term[0] = O*hf;

      compress(term[0],mpsxx::Left,cutoff[0]);
      compress(term[0],mpsxx::Right,cutoff[0]);

   }
   else{

      for(int i = term.size() - 1;i > 0;--i){

         term[i] = O * term[i - 1];

         compress(term[i],mpsxx::Left,cutoff[order[i] - 1]);
         compress(term[i],mpsxx::Right,cutoff[order[i] - 1]);

      }

      //remove the first element
      term.erase(term.begin());
      order.erase(order.begin());

   }

   ord++;

}

/**
 * evaluate the inproduct of a an operator with 2 eMPS obejcts, assume the operator connects only 2 particle terms
 */
double eMPS::inprod(const MPO<Quantum> &H,const MPS<Quantum> &hf,const eMPS &ket) const {

   if(ord == 0 && ket.ord == 0){//both are 0

      double E = mpsxx::inprod(mpsxx::Left,hf,H,hf) + mpsxx::inprod(mpsxx::Left,hf,H,term[0]) + mpsxx::inprod(mpsxx::Left,hf,H,ket.term[0]);

      for(int i = 0;i < cutoff.size();++i){

         E += mpsxx::inprod(mpsxx::Left,term[i],H,ket.term[i]);

         if(i + 1 < cutoff.size())
            E += mpsxx::inprod(mpsxx::Left,term[i],H,ket.term[i + 1]) + mpsxx::inprod(mpsxx::Left,ket.term[i],H,term[i + 1]);

      }

      return E;

   }
   else if(ord == 1 && ket.ord == 0){//one has no hf state

         double E = mpsxx::inprod(mpsxx::Left,hf,H,term[0]);

         for(int i = 0;i < cutoff.size();++i){

            E += mpsxx::inprod(mpsxx::Left,term[i],H,ket.term[i]);

            if(i + 1 < cutoff.size())
               E += mpsxx::inprod(mpsxx::Left,term[i],H,ket.term[i + 1]) + mpsxx::inprod(mpsxx::Left,ket.term[i],H,term[i + 1]);

         }

         return E;

   }
   else if(ord == 0 && ket.ord == 1){//ket has no hf

         double E = mpsxx::inprod(mpsxx::Left,hf,H,ket.term[0]);

         for(int i = 0;i < cutoff.size();++i){

            E += mpsxx::inprod(mpsxx::Left,term[i],H,ket.term[i]);

            if(i + 1 < cutoff.size())
               E += mpsxx::inprod(mpsxx::Left,term[i],H,ket.term[i + 1]) + mpsxx::inprod(mpsxx::Left,ket.term[i],H,term[i + 1]);

         }

         return E;

   }
   else{//the rest

      double E = 0.0;
      
      for(int i = 0;i < term.size();++i)//loop over bra
         for(int j = 0;j < ket.term.size();++j){//loop over ket

            if( abs(order[i] - ket.order[j]) < 2 )//if the difference in order is 0 or 1, Hamiltonian connects
               E += mpsxx::inprod(mpsxx::Left,term[i],H,ket.term[j]);

         }

      return E;

   }

}

/**
 * evaluate the dotproduct of two eMPS objects
 */
double eMPS::dot(const eMPS &ket) const {

   double N = 0.0;

   if(ord == 0 && ket.ord == 0)//both are 0
      N = 1.0;

   for(int i = 0;i < term.size();++i)//loop over bra
      for(int j = 0;j < ket.term.size();++j){//loop over ket

         if(order[i] == ket.order[j])
            N += mpsxx::dot(mpsxx::Left,term[i],ket.term[j]);

      }

   return N;

}
