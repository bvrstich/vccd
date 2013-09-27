#include <iostream>
#include <fstream>
#include <cmath>

using std::endl;
using std::ostream;
using std::ofstream;
using std::ifstream;

#include "include.h"

/**
 * standard constructor
 */
e_eMPS::e_eMPS(const MPO<Quantum> &O,const eMPS &ccd,const MPS<Quantum> &hf){

   term.resize(ccd.gcutoff().size());

   term[0] = ccd;
   term[0].mult(O,hf);

   for(int i = 1;i < term.size();++i){

      term[i] = term[i - 1];
      term[i].mult(O,hf);

   }

}

/**
 * copy constructor
 */
e_eMPS::e_eMPS(const e_eMPS &e_emps){

   term = e_emps.term;

}

/**
 * destructor
 */
e_eMPS::~e_eMPS(){ }

/**
 * fill the array E with the different inner products of the terms in this
 */
void e_eMPS::fillE(DArray<2> &Emat,const MPO<Quantum> &H,const MPS<Quantum> &hf,const eMPS &ccd) const{

   //first row
   for(int i = 0;i < term.size();++i)
      Emat(0,i + 1) = term[i].inprod(H,hf,ccd);

   for(int i = 0;i < term.size();++i)
      for(int j = i;j < term.size();++j)
         Emat(i + 1,j + 1) = term[i].inprod(H,hf,term[j]);

}

/**
 * fill the array E with the different inner products of the terms in this
 */
void e_eMPS::fillN(DArray<2> &Nmat,const eMPS &ccd) const{

   //first row
   for(int i = 0;i < term.size();++i)
      Nmat(0,i + 1) = term[i].dot(ccd);

   for(int i = 0;i < term.size();++i)
      for(int j = i;j < term.size();++j)
         Nmat(i + 1,j + 1) = term[i].dot(term[j]);

}
