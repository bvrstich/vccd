#include <iostream>
#include <iomanip>

using std::cout;
using std::endl;
using std::ostream;

#include "include.h"

namespace btas {

   /**
    * join two DArray objects together by blocking them together in a larger matrix C.
    * After this function is called C will contain:  [ A   0 ] 
    *                                                [ 0   B ]  
    * @param A input DArray<3> object
    * @param B input DArray<3> object
    * @param C output DArray<3> object, input will be destroyed by this function
    */
   void Djoin(const DArray<3> &A,const DArray<3> &B,DArray<3> &C){

      TinyVector<int,3> ashape = A.shape();
      TinyVector<int,3> bshape = B.shape();

      int d = ashape(1);

      TinyVector<int,3> cshape(ashape(0) + bshape(0),d,ashape(2) + bshape(2));
      C.resize(cshape);
      
      C = 0.0;

      for(int s = 0;s < d;++s){

         for(int i = 0;i < ashape(0);++i)
            for(int j = 0;j < ashape(2);++j)
               C(i,s,j) = A(i,s,j);

         for(int i = ashape(0);i < cshape(0);++i)
            for(int j = ashape(2);j < cshape(2);++j)
               C(i,s,j) = B(i - ashape(0),s,j - ashape(2));
      }

   }

   /**
    * join two QSDArray objects together by blocking them together in a larger matrix C.
    * After this function is called C will contain:  [ A   0 ] 
    *                                                [ 0   B ]  
    * @param A input QSDArray<3> object
    * @param B input QSDArray<3> object
    * @param C output QSDArray<3> object, input will be destroyed by this function
    */
   void QSDjoin(const QSDArray<3> &A,const QSDArray<3> &B,QSDArray<3> &C){

      //first calculate the new quantumnumbers and dimensions of C
      TinyVector<Qshapes,3> qa = A.qshape();
      TinyVector<Dshapes,3> da = A.dshape();
      TinyVector<Qshapes,3> qb = B.qshape();
      TinyVector<Dshapes,3> db = B.dshape();

      TinyVector<Qshapes,3> qc(qa);
      TinyVector<Dshapes,3> dc(da);

      //first 0 index
      Qshapes::iterator jt = qb[0].begin();

      Dshapes::iterator kt = dc[0].begin();
      Dshapes::iterator lt = db[0].begin();

      for(Qshapes::iterator it = qc[0].begin();it != qc[0].end();++it){

         while( (*jt) < (*it) && jt != qb[0].end() ){

            qc[0].insert(it,*jt);
            dc[0].insert(kt,*lt);

            ++jt;
            ++lt;

         }

         if( (*jt) == (*it) && jt != qb[0].end()  ){

            *kt += *lt;

            ++jt;
            ++lt;

         }

         ++kt;

      }

      //then 2 index
      jt = qb[2].begin();

      kt = dc[2].begin();
      lt = db[2].begin();

      for(Qshapes::iterator it = qc[2].begin();it != qc[2].end();++it){

         while( (*jt) < (*it) && jt != qb[2].end() ){

            qc[2].insert(it,*jt);
            dc[2].insert(kt,*lt);

            ++lt;
            ++jt;

         }

         if( (*jt) == (*it) && jt != qb[2].end()  ){

            *kt += *lt;

            ++jt;
            ++lt;

         }

         ++kt;

      }

      C.resize(Quantum::zero(),qc,dc);

   }

}
