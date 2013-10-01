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
e_eMPS::e_eMPS(const MPO<Quantum> &O,const eMPS &ccd,const MPS<Quantum> &hf) : std::vector< eMPS > (ccd.gcutoff().size()) {

   (*this)[0] = ccd;
   (*this)[0].mult_exc(O,hf);

   for(int i = 1;i < (*this).size();++i){

      (*this)[i] = (*this)[i - 1];
      (*this)[i].mult_exc(O,hf);

   }

}

/**
 * copy constructor
 */
e_eMPS::e_eMPS(const e_eMPS &e_emps) : std::vector< eMPS > (e_emps){ }

/**
 * destructor
 */
e_eMPS::~e_eMPS(){ }

/**
 * fill the array E with the different inner products of the (*this)s in this
 */
void e_eMPS::fillE(DArray<2> &Emat,const MPO<Quantum> &H,const MPS<Quantum> &hf,const eMPS &ccd) const{

   //first row
   for(int i = 0;i < (*this).size();++i)
      Emat(0,i + 1) = (*this)[i].inprod(H,hf,ccd);

   for(int i = 0;i < (*this).size();++i)
      for(int j = i;j < (*this).size();++j)
         Emat(i + 1,j + 1) = (*this)[i].inprod(H,hf,(*this)[j]);

}

/**
 * fill the array E with the different inner products of the (*this)s in this
 */
void e_eMPS::fillN(DArray<2> &Nmat,const eMPS &ccd) const{

   //first row
   for(int i = 0;i < (*this).size();++i)
      Nmat(0,i + 1) = (*this)[i].dot(ccd);

   for(int i = 0;i < (*this).size();++i)
      for(int j = i;j < (*this).size();++j)
         Nmat(i + 1,j + 1) = (*this)[i].dot((*this)[j]);

}
