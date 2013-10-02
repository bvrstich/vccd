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
eMPS::eMPS() : std::vector< MPS<Quantum> > (){ }

/**
 * standard constructor, initializes the T and hf
 */
eMPS::eMPS(const MPO<Quantum> &T,const MPS<Quantum> &hf,const std::vector<int> &cutoff_in) : std::vector< MPS<Quantum> >(cutoff_in.size()) {

   cutoff = cutoff_in;

   //form the list of contributing (*this)s in the expansion
   (*this)[0] = T*hf;

   compress((*this)[0],mpsxx::Left,cutoff[0]);
   compress((*this)[0],mpsxx::Right,cutoff[0]);

   for(int i = 1;i < cutoff.size();++i){

      (*this)[i] = T* (*this)[i - 1];
      compress((*this)[i],mpsxx::Left,cutoff[i]);
      compress((*this)[i],mpsxx::Right,cutoff[i]);
      scal(1.0/(i + 1.0),(*this)[i]);

   }

   order.resize(cutoff.size());

   for(int i = 0;i < order.size();++i)
      order[i] = i + 1;

   ord = 0;

}

/**
 * copy constructor
 */
eMPS::eMPS(const eMPS &emps) : std::vector< MPS<Quantum> >(emps) { 

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
 * evaluate the expectation value of an operator, assume the operator connects only 2 particle (*this)s
 */
double eMPS::eval(const MPO<Quantum> &H,const MPS<Quantum> &hf) const {

   double E = mpsxx::inprod(mpsxx::Left,hf,H,hf) + 2.0 * mpsxx::inprod(mpsxx::Left,hf,H,(*this)[0]);

   for(int i = 0;i < cutoff.size();++i){

      E += mpsxx::inprod(mpsxx::Left,(*this)[i],H,(*this)[i]);

      if(i + 1 < cutoff.size())
         E += 2.0 * mpsxx::inprod(mpsxx::Left,(*this)[i],H,(*this)[i + 1]);

   }

   return E;

}

/**
 * evaluate the expectation value of an operator, assume the operator connects only 2 particle (*this)s
 */
double eMPS::norm() const {

   double nrm = 1.0;

   for(int i = 0;i < cutoff.size();++i)
       nrm += mpsxx::dot(mpsxx::Left,(*this)[i],(*this)[i]);

   return nrm;

}

/**
 * @param hf hartree fock reference state
 * @param rank order of the expansion, if 1 MPS = hf + T hf , etc.: has to be larger than 0!
 * @param D dimension of the compression of the final MPS
 * @return the exponential as one MPS
 */
MPS<Quantum> eMPS::expand(const MPS<Quantum> &hf,int rank,int D) const{

   if(rank > 0){

   MPS<Quantum> sum = (*this)[0];

   //now sum all the terms together:
   for(int i = 1;i < rank;++i){

      axpy(1.0,(*this)[i],sum);
      compress(sum,mpsxx::Left,D);
      compress(sum,mpsxx::Right,D);

   }

   axpy(1.0,hf,sum);

   return sum;

   }
   else
      return hf;

}

/**
 * update the eMPS for a new MPO form of the T operator
 */
void eMPS::update(const MPO<Quantum> &T,const MPS<Quantum> &hf){

   //form the list of contributing (*this)s in the expansion
   (*this)[0] = T*hf;

   compress((*this)[0],mpsxx::Left,cutoff[0]);
   compress((*this)[0],mpsxx::Right,cutoff[0]);

   for(int i = 1;i < cutoff.size();++i){

      (*this)[i] = T*(*this)[i - 1];
      compress((*this)[i],mpsxx::Left,cutoff[i]);
      compress((*this)[i],mpsxx::Right,cutoff[i]);
      scal(1.0/(i + 1.0),(*this)[i]);

   }

}

/**
 * multiply this with an MPO form of a T2 excitation operator with different t amplitudes
 * @param O input MPO
 */
void eMPS::mult_exc(const MPO<Quantum> &O,const MPS<Quantum> &hf){

   if(ord == 0){

      for(int i = (*this).size() - 1;i > 0;--i){

         (*this)[i] = O * (*this)[i - 1];

         compress((*this)[i],mpsxx::Left,cutoff[i]);
         compress((*this)[i],mpsxx::Right,cutoff[i]);

      }

      (*this)[0] = O*hf;

      compress((*this)[0],mpsxx::Left,cutoff[0]);
      compress((*this)[0],mpsxx::Right,cutoff[0]);

   }
   else{

      for(int i = (*this).size() - 1;i > 0;--i){

         (*this)[i] = O * (*this)[i - 1];

         compress((*this)[i],mpsxx::Left,cutoff[order[i] - 1]);
         compress((*this)[i],mpsxx::Right,cutoff[order[i] - 1]);

      }

      //remove the first element
      (*this).erase((*this).begin());
      order.erase(order.begin());

   }

   ord++;

}

/**
 * evaluate the inproduct of a an operator with 2 eMPS obejcts, assume the operator connects only 2 particle (*this)s
 */
double eMPS::inprod(const MPO<Quantum> &H,const MPS<Quantum> &hf,const eMPS &ket) const {

   if(ord == 0 && ket.ord == 0){//both are 0

      double E = mpsxx::inprod(mpsxx::Left,hf,H,hf) + mpsxx::inprod(mpsxx::Left,hf,H,(*this)[0]) + mpsxx::inprod(mpsxx::Left,hf,H,ket[0]);

      for(int i = 0;i < cutoff.size();++i){

         E += mpsxx::inprod(mpsxx::Left,(*this)[i],H,ket[i]);

         if(i + 1 < cutoff.size())
            E += mpsxx::inprod(mpsxx::Left,(*this)[i],H,ket[i + 1]) + mpsxx::inprod(mpsxx::Left,ket[i],H,(*this)[i + 1]);

      }

      return E;

   }
   else if(ord == 1 && ket.ord == 0){//one has no hf state

      double E = mpsxx::inprod(mpsxx::Left,hf,H,(*this)[0]);

      for(int i = 0;i < cutoff.size();++i){

         E += mpsxx::inprod(mpsxx::Left,(*this)[i],H,ket[i]);

         if(i + 1 < cutoff.size())
            E += mpsxx::inprod(mpsxx::Left,(*this)[i],H,ket[i + 1]) + mpsxx::inprod(mpsxx::Left,ket[i],H,(*this)[i + 1]);

      }

      return E;

   }
   else if(ord == 0 && ket.ord == 1){//ket has no hf

      double E = mpsxx::inprod(mpsxx::Left,hf,H,ket[0]);

      for(int i = 0;i < cutoff.size();++i){

         E += mpsxx::inprod(mpsxx::Left,(*this)[i],H,ket[i]);

         if(i + 1 < cutoff.size())
            E += mpsxx::inprod(mpsxx::Left,(*this)[i],H,ket[i + 1]) + mpsxx::inprod(mpsxx::Left,ket[i],H,(*this)[i + 1]);

      }

      return E;

   }
   else{//the rest

      double E = 0.0;

      for(int i = 0;i < (*this).size();++i)//loop over bra
         for(int j = 0;j < ket.size();++j){//loop over ket

            if( abs(order[i] - ket.order[j]) < 2 )//if the difference in order is 0 or 1, Hamiltonian connects
               E += mpsxx::inprod(mpsxx::Left,(*this)[i],H,ket[j]);

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

   for(int i = 0;i < (*this).size();++i)//loop over bra
      for(int j = 0;j < ket.size();++j){//loop over ket

         if(order[i] == ket.order[j])
            N += mpsxx::dot(mpsxx::Left,(*this)[i],ket[j]);

      }

   return N;

}

/**
 * print the dimensions of the mps of order order_in
 * @param order_in input order
 */
void eMPS::print_dim(int order_in){

   for(int i = 0;i < (*this)[order_in].size();++i){

      cout << endl;
      cout << i << endl;
      cout << endl;
      cout << (*this)[order_in][i].qshape() << endl;
      cout << (*this)[order_in][i].dshape() << endl;
      cout << endl;

   }

}

/**
 * print the dimensions of the mps of order order_in
 * @param order_in input order
 */
void eMPS::print_tot_dim(int order_in){

   for(int i = 0;i < (*this)[order_in].size();++i){

      cout << i << "\t";

      //row dim
      int tmp = 0;

      for(int j = 0;j < (*this)[order_in][i].dshape(0).size();++j)
         tmp += (*this)[order_in][i].dshape(0)[j];

      cout << tmp << "\t";

      //col dim
      tmp = 0;

      for(int j = 0;j < (*this)[order_in][i].dshape(2).size();++j)
         tmp += (*this)[order_in][i].dshape(2)[j];

      cout << tmp << endl;

   }

}
