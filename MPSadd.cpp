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
    * will copy a DArray object in a larger matrix C.
    * After this function is called C will contain:  [ A   0 ] 
    *                                                [ 0   0 ]  
    * @param A input DArray<3> object
    * @param C output DArray<3> object, input will be destroyed by this function
    */
   void Djoin_upper(const DArray<3> &A,DArray<3> &C){

      TinyVector<int,3> ashape = A.shape();

      int d = ashape(1);

      TinyVector<int,3> cshape = C.shape();
      
      C = 0.0;

      for(int s = 0;s < d;++s)
         for(int i = 0;i < ashape(0);++i)
            for(int j = 0;j < ashape(2);++j)
               C(i,s,j) = A(i,s,j);

   }

   /**
    * will copy a DArray object in a larger matrix C.
    * After this function is called C will contain:  [ 0   0 ] 
    *                                                [ 0   B ]  
    * @param B input DArray<3> object
    * @param C output DArray<3> object, input will be destroyed by this function
    */
   void Djoin_lower(const DArray<3> &B,DArray<3> &C){

      TinyVector<int,3> bshape = B.shape();

      int d = bshape(1);

      TinyVector<int,3> cshape = C.shape();
      
      C = 0.0;

      for(int s = 0;s < d;++s)
         for(int i = 0;i < bshape(0);++i)
            for(int j = 0;j < bshape(2);++j)
               C(cshape(0) - bshape(0) + i,s,cshape(2) - bshape(2) + j) = B(i,s,j);

   }

   /**
    * join two 'left boundary' DArray objects together by blocking them together in a larger matrix C.
    * Only the right index will increase in size!
    * After this function is called C will contain:  [ A  B ]
    * @param A input DArray<3> object
    * @param B input DArray<3> object
    * @param C output DArray<3> object, input will be destroyed by this function
    */
   void Djoin_ledge(const DArray<3> &A,const DArray<3> &B,DArray<3> &C){

      TinyVector<int,3> ashape = A.shape();
      TinyVector<int,3> bshape = B.shape();

      if(ashape(0) != 1)
         cout << "No left edge!" << endl;

      if(bshape(0) != 1)
         cout << "No left edge!" << endl;

      int d = ashape(1);

      TinyVector<int,3> cshape(1,d,ashape(2) + bshape(2));
      C.resize(cshape);

      C = 0.0;

      for(int s = 0;s < d;++s){

         for(int j = 0;j < ashape(2);++j)
            C(0,s,j) = A(0,s,j);

         for(int j = ashape(2);j < cshape(2);++j)
            C(0,s,j) = B(0,s,j - ashape(2));

      }

   }

   /**
    * join two 'right boundary' DArray objects together by blocking them together in a larger matrix C.
    * Only the left index will increase in size!
    * After this function is called C will contain:  [ A ]
    *                                                [ B ]
    * @param A input DArray<3> object
    * @param B input DArray<3> object
    * @param C output DArray<3> object, input will be destroyed by this function
    */
   void Djoin_redge(const DArray<3> &A,const DArray<3> &B,DArray<3> &C){

      TinyVector<int,3> ashape = A.shape();
      TinyVector<int,3> bshape = B.shape();

      if(ashape(2) != 1)
         cout << "No right edge!" << endl;

      if(bshape(2) != 1)
         cout << "No right edge!" << endl;

      int d = ashape(1);

      TinyVector<int,3> cshape(ashape(0) + bshape(0),d,1);
      C.resize(cshape);

      C = 0.0;

      for(int s = 0;s < d;++s){

         for(int i = 0;i < ashape(0);++i)
            C(i,s,0) = A(i,s,0);

         for(int i = ashape(0);i < cshape(0);++i)
            C(i,s,0) = B(i - ashape(0),s,0);

      }

   }

   /**
    * join two 'left boundary' DArray objects together by blocking them together in a larger matrix C.
    * Only the right index will increase in size! no resizing the C array
    * After this function is called C will contain:  [ A  B ]
    * @param A input DArray<3> object
    * @param B input DArray<3> object
    * @param C output DArray<3> object, input will be destroyed by this function, assume dimensions are correct on input
    */
   void Djoin_ledge_no_reshape(const DArray<3> &A,const DArray<3> &B,DArray<3> &C){

      TinyVector<int,3> ashape = A.shape();
      TinyVector<int,3> bshape = B.shape();
      TinyVector<int,3> cshape = C.shape();

      if(ashape(0) != 1)
         cout << "No left edge!" << endl;

      if(bshape(0) != 1)
         cout << "No left edge!" << endl;

      int d = ashape(1);

      C = 0.0;

      for(int s = 0;s < d;++s){

         for(int j = 0;j < ashape(2);++j)
            C(0,s,j) = A(0,s,j);

         for(int j = ashape(2);j < cshape(2);++j)
            C(0,s,j) = B(0,s,j - ashape(2));

      }

   }

   /**
    * join two 'right boundary' DArray objects together by blocking them together in a larger matrix C.
    * Only the left index will increase in size! No resizing the C array
    * After this function is called C will contain:  [ A ]
    *                                                [ B ]
    * @param A input DArray<3> object
    * @param B input DArray<3> object
    * @param C output DArray<3> object, input will be destroyed by this function, assume dimensions are correct on input
    */
   void Djoin_redge_no_reshape(const DArray<3> &A,const DArray<3> &B,DArray<3> &C){

      TinyVector<int,3> ashape = A.shape();
      TinyVector<int,3> bshape = B.shape();
      TinyVector<int,3> cshape = C.shape();

      if(ashape(2) != 1)
         cout << "No right edge!" << endl;

      if(bshape(2) != 1)
         cout << "No right edge!" << endl;

      int d = ashape(1);

      C = 0.0;

      for(int s = 0;s < d;++s){

         for(int i = 0;i < ashape(0);++i)
            C(i,s,0) = A(i,s,0);

         for(int i = ashape(0);i < cshape(0);++i)
            C(i,s,0) = B(i - ashape(0),s,0);

      }

   }

   /**
    * join two DArray objects together by blocking them together in a larger matrix C. No reshaping of matrix C!
    * After this function is called C will contain:  [ A   0 ] 
    *                                                [ 0   B ]  
    * @param A input DArray<3> object
    * @param B input DArray<3> object
    * @param C output DArray<3> object, input will be destroyed by this function. function assumes C already has the right dimensions!
    */
   void Djoin_no_reshape(const DArray<3> &A,const DArray<3> &B,DArray<3> &C){

      TinyVector<int,3> ashape = A.shape();
      TinyVector<int,3> bshape = B.shape();
      TinyVector<int,3> cshape = C.shape();

      int d = ashape(1);

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
      TinyVector<Qshapes,3> qb = B.qshape();

      TinyVector<Qshapes,3> qc(qa);

      TinyVector<Dshapes,3> da = A.dshape();
      TinyVector<Dshapes,3> db = B.dshape();

      TinyVector<Dshapes,3> dc(da);

      //first 0 index
      Qshapes::iterator jt = qb[0].begin();
      Dshapes::iterator lt = db[0].begin();

      Dshapes::iterator kt = dc[0].begin();

      int i = 0;

      for(Qshapes::iterator it = qc[0].begin();it != qc[0].end();++it){

         while( (*jt) < (*it) ){

            qc[0].insert(it,*jt);
            dc[0].insert(kt,*lt);

            //ridiculous but I need a new iterator
            ++i;

            it = qc[0].begin() + i;
            kt = dc[0].begin() + i;

            if(it == qc[0].end())
               break;

            ++jt;
            ++lt;

            if(jt == qb[0].end())
               break;

         }

         if(it != qc[0].end() && jt != qb[0].end() ){

            if( (*jt) == (*it) ){

               *kt += *lt;

               ++jt;
               ++lt;

            }

         }

         if(jt == qb[0].end())
            break;

         ++i;
         ++kt;

      }

      //now add stuff to the end
      while(jt != qb[0].end()){

         qc[0].push_back(*jt);
         dc[0].push_back(*lt);

         ++jt;
         ++lt;

      }

      //then 2 index: this runs from high to low
      jt = qb[2].begin();
      lt = db[2].begin();

      kt = dc[2].begin();

      i = 0;

      for(Qshapes::iterator it = qc[2].begin();it != qc[2].end();++it){

         while( (*jt) > (*it) ){

            qc[2].insert(it,*jt);
            dc[2].insert(kt,*lt);

            //ridiculous but I need a new iterator
            ++i;

            it = qc[2].begin() + i;
            kt = dc[2].begin() + i;

            if(it == qc[2].end())
               break;

            ++jt;
            ++lt;

            if(jt == qb[2].end())
               break;

         }

         if( it != qc[2].end() && jt != qb[2].end() ){

            if( (*jt) == (*it) ){

               *kt += *lt;

               ++lt;
               ++jt;

            }

         }

         if(jt == qb[2].end())
            break;

         ++i;
         ++kt;

      }

      //now add stuff to the end
      while(jt != qb[2].end()){

         qc[2].push_back(*jt);
         dc[2].push_back(*lt);

         ++jt;
         ++lt;

      }

      C.resize(Quantum::zero(),qc,dc);

      Qshapes qindc(3);

      for(SDArray<3>::const_iterator itc = C.begin();itc != C.end();++itc){

         //put the non-zero index to quantum numbers
         qindex(C,itc->first,qindc);

         //check if the same numbers are present in A
         Qshapes qinda(3);

         int flag_A = 0;

         SDArray<3>::const_iterator A_loc;

         for(SDArray<3>::const_iterator ita = A.begin();ita != A.end();++ita){

            qindex(A,ita->first,qinda);

            if(qindc == qinda){

               A_loc = ita;
               flag_A = 1;

            }

         }

         //check if the same numbers are present in B
         Qshapes qindb(3);

         int flag_B = 0;

         SDArray<3>::const_iterator B_loc;

         for(SDArray<3>::const_iterator itb = B.begin();itb != B.end();++itb){

            qindex(B,itb->first,qindb);

            if(qindc == qindb){

               B_loc = itb;
               flag_B = 1;

            }

         }

         if(flag_A == 1){//if quantumnumber is present in A

            if(flag_B == 1)//and quantumnumber is present in B
               Djoin_no_reshape(*(A_loc->second),*(B_loc->second),*(itc->second));
            else//copy A to C
               Djoin_upper(*(A_loc->second),*(itc->second));

         }
         else{//if quantumnumber is not present in A

            if(flag_B == 1)//but it is in B: copy B to C
               Djoin_lower(*(B_loc->second),*(itc->second));
            else//this block is not present in both A and B: set to zero
               *(itc->second) = 0.0;

         }

      }

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

      //first calculate the new quantumnumbers and dimensions of C
      TinyVector<Qshapes,3> qa = A.qshape();
      TinyVector<Qshapes,3> qb = B.qshape();

      TinyVector<Qshapes,3> qc(qa);

      TinyVector<Dshapes,3> da = A.dshape();
      TinyVector<Dshapes,3> db = B.dshape();

      TinyVector<Dshapes,3> dc(da);

      //only 2 index: this runs from high to low
      Qshapes::const_iterator jt = qb[2].begin();
      Dshapes::const_iterator lt = db[2].begin();

      Dshapes::iterator kt = dc[2].begin();

      int i = 0;

      for(Qshapes::iterator it = qc[2].begin();it != qc[2].end();++it){

         while( (*jt) > (*it) ){

            qc[2].insert(it,*jt);
            dc[2].insert(kt,*lt);

            //ridiculous but I need a new iterator
            ++i;

            it = qc[2].begin() + i;
            kt = dc[2].begin() + i;

            if(it == qc[2].end())
               break;

            ++jt;

            if(jt == qb[2].end())
               break;

         }

         if( it != qc[2].end() && jt != qb[2].end() ){

            if( (*jt) == (*it) ){

               *kt += *lt;

               ++lt;
               ++jt;

            }

         }

         if(jt == qb[2].end())
            break;

         ++i;
         ++kt;

      }

      //now add stuff to the end
      while(jt != qb[2].end()){

         qc[2].push_back(*jt);
         dc[2].push_back(*lt);

         ++jt;
         ++lt;

      }

      C.resize(Quantum::zero(),qc,dc);

      Qshapes qindc(3);

      for(SDArray<3>::const_iterator itc = C.begin();itc != C.end();++itc){

         //put the non-zero index to quantum numbers
         qindex(C,itc->first,qindc);

         //check if the same numbers are present in A
         Qshapes qinda(3);

         int flag_A = 0;

         SDArray<3>::const_iterator A_loc;

         for(SDArray<3>::const_iterator ita = A.begin();ita != A.end();++ita){

            qindex(A,ita->first,qinda);

            if(qindc == qinda){

               A_loc = ita;
               flag_A = 1;

            }

         }

         //check if the same numbers are present in B
         Qshapes qindb(3);

         int flag_B = 0;

         SDArray<3>::const_iterator B_loc;

         for(SDArray<3>::const_iterator itb = B.begin();itb != B.end();++itb){

            qindex(B,itb->first,qindb);

            if(qindc == qindb){

               B_loc = itb;
               flag_B = 1;

            }

         }

         if(flag_A == 1){//if quantumnumber is present in A

            if(flag_B == 1)//and quantumnumber is present in B
               Djoin_ledge_no_reshape(*(A_loc->second),*(B_loc->second),*(itc->second));
            else//copy A to C
               Dcopy(*(A_loc->second),*(itc->second));

         }
         else{//if quantumnumber is not present in A

            if(flag_B == 1)//but it is in B: copy B to C
               Dcopy(*(B_loc->second),*(itc->second));
            else//this block is not present in both A and B: set to zero
               *(itc->second) = 0.0;

         }

      }

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

      //first calculate the new quantumnumbers and dimensions of C
      TinyVector<Qshapes,3> qa = A.qshape();
      TinyVector<Qshapes,3> qb = B.qshape();

      TinyVector<Qshapes,3> qc(qa);

      TinyVector<Dshapes,3> da = A.dshape();
      TinyVector<Dshapes,3> db = B.dshape();

      TinyVector<Dshapes,3> dc(da);

      //first 0 index
      Qshapes::iterator jt = qb[0].begin();
      Dshapes::iterator lt = db[0].begin();

      Dshapes::iterator kt = dc[0].begin();

      int i = 0;

      for(Qshapes::iterator it = qc[0].begin();it != qc[0].end();++it){

         while( (*jt) < (*it) ){

            qc[0].insert(it,*jt);
            dc[0].insert(kt,*lt);

            //ridiculous but I need a new iterator
            ++i;

            it = qc[0].begin() + i;
            kt = dc[0].begin() + i;

            if(it == qc[0].end())
               break;

            ++jt;
            ++lt;

            if(jt == qb[0].end())
               break;

         }

         if(it != qc[0].end() && jt != qb[0].end() ){

            if( (*jt) == (*it) ){

               *kt += *lt;

               ++jt;
               ++lt;

            }

         }

         if(jt == qb[0].end())
            break;

         ++i;
         ++kt;

      }

      //now add stuff to the end
      while(jt != qb[0].end()){

         qc[0].push_back(*jt);
         dc[0].push_back(*lt);

         ++jt;
         ++lt;

      }

      C.resize(Quantum::zero(),qc,dc);

      Qshapes qindc(3);

      for(SDArray<3>::const_iterator itc = C.begin();itc != C.end();++itc){

         //put the non-zero index to quantum numbers
         qindex(C,itc->first,qindc);

         //check if the same numbers are present in A
         Qshapes qinda(3);

         int flag_A = 0;

         SDArray<3>::const_iterator A_loc;

         for(SDArray<3>::const_iterator ita = A.begin();ita != A.end();++ita){

            qindex(A,ita->first,qinda);

            if(qindc == qinda){

               A_loc = ita;
               flag_A = 1;

            }

         }

         //check if the same numbers are present in B
         Qshapes qindb(3);

         int flag_B = 0;

         SDArray<3>::const_iterator B_loc;

         for(SDArray<3>::const_iterator itb = B.begin();itb != B.end();++itb){

            qindex(B,itb->first,qindb);

            if(qindc == qindb){

               B_loc = itb;
               flag_B = 1;

            }

         }

         if(flag_A == 1){//if quantumnumber is present in A

            if(flag_B == 1)//and quantumnumber is present in B
               Djoin_redge_no_reshape(*(A_loc->second),*(B_loc->second),*(itc->second));
            else//copy A to upper block of C
               Djoin_upper(*(A_loc->second),*(itc->second));

         }
         else{//if quantumnumber is not present in A

            if(flag_B == 1)//but it is in B: copy B to lower block of C
               Djoin_lower(*(B_loc->second),*(itc->second));
            else//this block is not present in both A and B: set to zero
               *(itc->second) = 0.0;

         }

      }

   }


   /** 
    * given a block_tag from the SDArray::iterator
    * @param A input QSDArray object
    * @param block_tag index of the non-zero block (key of the map)
    * @param q Qshapes of length 3 will contain the 3 indices at exit
    * @return the Quantumnumbers belonging to the non-zero block 
    */
   void qindex(const QSDArray<3> &A,int block_tag,Qshapes &q){

      TinyVector<int,3> ind = A.index(block_tag);

      for(int i = 0;i < 3;++i)
         q[i] = A.qshape(i)[ind[i]];

   }

}
