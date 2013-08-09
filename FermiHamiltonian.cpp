#include <iostream>
#include <iomanip>

using std::cout;
using std::endl;
using std::ostream;

#include "include.h"

namespace btas {
   
   /**
    * construct an MPO which creates a fermion at site 'site'
    */
   MPO creator(int L,int d,int site){

      MPO mpo(L);

      //first set the quantumnumbers, before
      Qshapes<Quantum> qp;
      physical(d,qp);

      Qshapes<Quantum> qz; // 0 quantum number
      qz.push_back(Quantum(0));

      Qshapes<Quantum> qt; // total quantum number
      qt.push_back(Quantum(1));

      //resize & set to 0
      for(int i = 0; i < site; ++i)
         mpo[i].resize(Quantum::zero(),make_array(qz,qp,-qp,qz));

      //the quantumnumbers on the site of the operation
      mpo[site].resize(Quantum::zero(),make_array(qz,qp,-qp,-qt));

      //the quantumnumbers after
      for(int i = site + 1;i < L;++i)
         mpo[i].resize(Quantum::zero(),make_array(qt,qp,-qp,-qt));

      DArray<4> Ip(1,1,1,1);
      Ip = 1;

      DArray<4> Im(1,1,1,1);
      Im = -1;

      //fill it up before
      for(int i = 0;i < site;++i){

         mpo[i].insert(shape(0,0,0,0),Ip);
         mpo[i].insert(shape(0,1,1,0),Im);

      }

      //on
      mpo[site].insert(shape(0,1,0,0),Ip);

      //after the operator
      for(int i = site + 1;i < L;++i){

         mpo[i].insert(shape(0,0,0,0),Ip);
         mpo[i].insert(shape(0,1,1,0),Ip);

      }

      return mpo;

   }

   /**
    * construct an MPO which annihilates a fermion at site 'site'
    */
   MPO annihilator(int L,int d,int site){

      MPO mpo(L);

      //first set the quantumnumbers, before
      Qshapes<Quantum> qp;
      physical(d,qp);

      Qshapes<Quantum> qz; // 0 quantum number
      qz.push_back(Quantum(0));

      Qshapes<Quantum> qt; // total quantum number
      qt.push_back(Quantum(-1));

      //resize & set to 0
      for(int i = 0; i < site; ++i)
         mpo[i].resize(Quantum::zero(),make_array(qz,qp,-qp,qz));

      //the quantumnumbers on the site of the operation
      mpo[site].resize(Quantum::zero(),make_array(qz,qp,-qp,-qt));

      //the quantumnumbers after
      for(int i = site + 1;i < L;++i)
         mpo[i].resize(Quantum::zero(),make_array(qt,qp,-qp,-qt));

      DArray<4> Ip(1,1,1,1);
      Ip = 1;

      DArray<4> Im(1,1,1,1);
      Im = -1;

      //fill it up before
      for(int i = 0;i < site;++i){

         mpo[i].insert(shape(0,0,0,0),Ip);
         mpo[i].insert(shape(0,1,1,0),Im);

      }

      //on
      mpo[site].insert(shape(0,0,1,0),Ip);

      //after the operator
      for(int i = site + 1;i < L;++i){

         mpo[i].insert(shape(0,0,0,0),Ip);
         mpo[i].insert(shape(0,1,1,0),Ip);

      }

      return mpo;

   }

}
