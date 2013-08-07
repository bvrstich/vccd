#include <iostream>
#include <iomanip>

using std::cout;
using std::endl;
using std::ostream;

#include "include.h"

namespace btas {

   /**
    * join two QSDArray objects together by blocking them together in a larger matrix C.
    * After this function is called C will contain:  [ A   0 ] 
    *                                                [ 0   B ]  
    * @param A input QSDArray<3> object
    * @param B input QSDArray<3> object
    * @param C output QSDArray<3> object, input will be destroyed by this function
    */
   void QSDjoin(const QSDArray<3> &A,const QSDArray<3> &B,QSDArray<3> &C){

      //first calculate the new quantumnumbers of C
      TVector<Qshapes<Quantum>,3> qa = A.qshape();
      TVector<Qshapes<Quantum>,3> qb = B.qshape();

      TVector<Qshapes<Quantum>,3> qc(qa);

      for(int i = 0;i < qb[0].size();++i)
         qc[0].push_back(qb[0][i]);

      for(int i = 0;i < qb[2].size();++i)
         qc[2].push_back(qb[2][i]);

      C.resize(Quantum::zero(),qc);

      //now insert the blocks
      int ald = qa[0].size();
      int ard = qa[2].size();

      int d = qa[1].size();

      //add A
      IVector<3> block_index;

      for(QSDArray<3>::const_iterator it = A.begin();it != A.end();++it){

         block_index = A.index(it->first);
         C.insert(block_index,*it->second);

      }

      //add B
      for(QSDArray<3>::const_iterator it = B.begin();it != B.end();++it){

         block_index = B.index(it->first);

         block_index[0] += ald;
         block_index[2] += ard;

         C.insert(block_index,*it->second);

      }

      //merge the row quantumnumbers together
      QSDArray<3> tmp;

      TVector<Qshapes<Quantum>,1> qmerge;
      TVector<Dshapes,1> dmerge;

      qmerge[0] = C.qshape(0);
      dmerge[0] = C.dshape(0);

      QSTmergeInfo<1> info(qmerge,dmerge);

      //then merge
      QSTmerge(info,C,tmp);

      //column quantumnumbers
      qmerge[0] = tmp.qshape(2);
      dmerge[0] = tmp.dshape(2);

      info.reset(qmerge,dmerge);

      //then merge
      QSTmerge(tmp,info,C);

   }

   /**
    * join two 'left boundary' QSDArray objects together by blocking them together in a larger matrix C.
    * Only the right index will change
    * After this function is called C will contain:  [ A  B ]  
    * @param A input QSDArray<3> object
    * @param B input QSDArray<3> object
    * @param C output QSDArray<3> object, input will be destroyed by this function
    */
   void QSDjoin_ledge(const QSDArray<3> &A,const QSDArray<3> &B,QSDArray<3> &C){

      //first calculate the new quantumnumbers of C
      TVector<Qshapes<Quantum>,3> qa = A.qshape();
      TVector<Qshapes<Quantum>,3> qb = B.qshape();

      TVector<Qshapes<Quantum>,3> qc(qa);

      for(int i = 0;i < qb[2].size();++i)
         qc[2].push_back(qb[2][i]);

      QSDArray<3> tmp;

      tmp.resize(Quantum::zero(),qc);

      //now insert the blocks
      int ard = qa[2].size();

      int d = qa[1].size();

      //add A
      IVector<3> block_index;

      for(QSDArray<3>::const_iterator it = A.begin();it != A.end();++it){

         block_index = A.index(it->first);
         tmp.insert(block_index,*it->second);

      }

      //add B
      for(QSDArray<3>::const_iterator it = B.begin();it != B.end();++it){

         block_index = B.index(it->first);

         block_index[2] += ard;

         tmp.insert(block_index,*it->second);

      }

      //merge the column quantumnumbers together
      TVector<Qshapes<Quantum>,1> qmerge;
      TVector<Dshapes,1> dmerge;

      qmerge[0] = tmp.qshape(2);
      dmerge[0] = tmp.dshape(2);

      QSTmergeInfo<1> info(qmerge,dmerge);

      //then merge
      QSTmerge(tmp,info,C);

   }

   /**
    * join two 'right boundary' QSDArray objects together by blocking them together in a larger matrix C.
    * After this function is called C will contain:  [ A ] 
    *                                                [ B ]  only 0 index will increase!
    * @param A input QSDArray<3> object
    * @param B input QSDArray<3> object
    * @param C output QSDArray<3> object, input will be destroyed by this function
    */
   void QSDjoin_redge(const QSDArray<3> &A,const QSDArray<3> &B,QSDArray<3> &C){

      //first calculate the new quantumnumbers of C
      TVector<Qshapes<Quantum>,3> qa = A.qshape();
      TVector<Qshapes<Quantum>,3> qb = B.qshape();

      TVector<Qshapes<Quantum>,3> qc(qa);

      for(int i = 0;i < qb[0].size();++i)
         qc[0].push_back(qb[0][i]);

      QSDArray<3> tmp;
      tmp.resize(Quantum::zero(),qc);

      //now insert the blocks
      int ald = qa[0].size();

      int d = qa[1].size();

      //add A
      IVector<3> block_index;

      for(QSDArray<3>::const_iterator it = A.begin();it != A.end();++it){

         block_index = A.index(it->first);
         tmp.insert(block_index,*it->second);

      }

      //add B
      for(QSDArray<3>::const_iterator it = B.begin();it != B.end();++it){

         block_index = B.index(it->first);

         block_index[0] += ald;

         tmp.insert(block_index,*it->second);

      }

      //merge the row quantumnumbers together
      TVector<Qshapes<Quantum>,1> qmerge;
      TVector<Dshapes,1> dmerge;

      qmerge[0] = tmp.qshape(0);
      dmerge[0] = tmp.dshape(0);

      QSTmergeInfo<1> info(qmerge,dmerge);

      //then merge into C
      QSTmerge(info,tmp,C);

   }

}
