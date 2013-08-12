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

   /**
    * construct an MPO which returns the local particle number operator a^+_i a_i
    */
   MPO n_loc(int L,int d,int site){

      MPO mpo(L);

      //first set the quantumnumbers, before
      Qshapes<Quantum> qp;
      physical(d,qp);

      Qshapes<Quantum> qz; // 0 quantum number
      qz.push_back(Quantum(0));

      //resize & set to 0
      for(int i = 0; i < L; ++i)
         mpo[i].resize(Quantum::zero(),make_array(qz,qp,-qp,qz));

      DArray<4> I(1,1,1,1);
      I = 1;

      //fill it up before
      for(int i = 0;i < site;++i){

         mpo[i].insert(shape(0,0,0,0),I);
         mpo[i].insert(shape(0,1,1,0),I);

      }

      //on
      mpo[site].insert(shape(0,1,1,0),I);

      //after the operator
      for(int i = site + 1;i < L;++i){

         mpo[i].insert(shape(0,0,0,0),I);
         mpo[i].insert(shape(0,1,1,0),I);

      }

      return mpo;

   }

   /**
    * @return MPO which returns the total particle number operator \sum_i a^+_i a_i
    */
   MPO N_tot(int L,int d){

      MPO mpo(L);

      //first set the quantumnumbers, before
      Qshapes<Quantum> qp;
      physical(d,qp);

      Qshapes<Quantum> qz; // 0 quantum number
      qz.push_back(Quantum(0));

      Qshapes<Quantum> qi;
      qi.push_back(Quantum::zero());
      qi.push_back(Quantum::zero());
      
      Qshapes<Quantum> qo;
      qo.push_back(Quantum::zero());
      qo.push_back(Quantum::zero());

      mpo[0].resize(Quantum::zero(),make_array(qz,qp,-qp,qo));

      //resize & set to 0
      for(int i = 1; i < L-1; ++i)
         mpo[i].resize(Quantum::zero(),make_array(qi,qp,-qp,qo));

      mpo[L-1].resize(Quantum::zero(),make_array(qi,qp,-qp,qz));

      DArray<4> I(1,1,1,1);
      I = 1;

      //left
      mpo[0].insert(shape(0,0,0,0),I);
      mpo[0].insert(shape(0,1,1,0),I);

      mpo[0].insert(shape(0,1,1,1),I);

      //middle
      for(int i = 1;i < L-1;++i){

         mpo[i].insert(shape(0,0,0,0),I);
         mpo[i].insert(shape(0,1,1,0),I);

         mpo[i].insert(shape(0,1,1,1),I);

         mpo[i].insert(shape(1,0,0,1),I);
         mpo[i].insert(shape(1,1,1,1),I);

      }

      //right
      mpo[L-1].insert(shape(0,1,1,0),I);

      mpo[L-1].insert(shape(1,0,0,0),I);
      mpo[L-1].insert(shape(1,1,1,0),I);

      //merge everything together
      TVector<Qshapes<Quantum>,1> qmerge;
      TVector<Dshapes,1> dmerge;

      qmerge[0] = mpo[0].qshape(3);
      dmerge[0] = mpo[0].dshape(3);

      QSTmergeInfo<1> info(qmerge,dmerge);

      QSDArray<4> tmp;
      QSTmerge(mpo[0],info,tmp);

      mpo[0] = tmp;

      for(int i = 1;i < L - 1;++i){

         //first merge the row
         qmerge[0] = mpo[i].qshape(0);
         dmerge[0] = mpo[i].dshape(0);

         info.reset(qmerge,dmerge);

         tmp.clear();

         QSTmerge(info,mpo[i],tmp);

         //then merge the column
         qmerge[0] = tmp.qshape(3);
         dmerge[0] = tmp.dshape(3);

         info.reset(qmerge,dmerge);

         mpo[i].clear();

         QSTmerge(tmp,info,mpo[i]);

      }

      //only merge row for i = L - 1
      qmerge[0] = mpo[L - 1].qshape(0);
      dmerge[0] = mpo[L - 1].dshape(0);

      info.reset(qmerge,dmerge);

      tmp.clear();

      QSTmerge(info,mpo[L - 1],tmp);

      mpo[L - 1] = tmp;

      return mpo;

   }


   /**
    * @return MPO which contains the nearest neighbour hopping operator
    */
   MPO hopping(int L,int d){

      MPO mpo(L);

      //first set the quantumnumbers, before
      Qshapes<Quantum> qp;
      physical(d,qp);

      Qshapes<Quantum> qz; // 0 quantum number
      qz.push_back(Quantum(0));

      Qshapes<Quantum> qi;
      qi.push_back(Quantum::zero());//0
      qi.push_back(Quantum(1));//a
      qi.push_back(Quantum(-1));//a^+
      qi.push_back(Quantum::zero());//I
      
      Qshapes<Quantum> qo;
      qo.push_back(Quantum::zero());//I
      qo.push_back(Quantum(-1));//a^+
      qo.push_back(Quantum(1));//a
      qo.push_back(Quantum::zero());//0

      mpo[0].resize(Quantum::zero(),make_array(qz,qp,-qp,qo));

      //resize & set to 0
      for(int i = 1; i < L-1; ++i)
         mpo[i].resize(Quantum::zero(),make_array(qi,qp,-qp,qo));

      mpo[L-1].resize(Quantum::zero(),make_array(qi,qp,-qp,qz));

      DArray<4> I(1,1,1,1);
      I = 1;

      DArray<4> O(1,1,1,1);
      O = 0;

      //left
      mpo[0].insert(shape(0,0,0,0),I);
      mpo[0].insert(shape(0,1,1,0),I);

      mpo[0].insert(shape(0,1,0,1),I);//a^+

      mpo[0].insert(shape(0,0,1,2),I);//a

      mpo[0].insert(shape(0,0,0,3),O);//zero, you have to put something here: else merge will reduce the rank!

      //middle
      for(int i = 1;i < L-1;++i){

         mpo[i].insert(shape(0,0,0,0),I);
         mpo[i].insert(shape(0,1,1,0),I);

         mpo[i].insert(shape(0,1,0,1),I);//a^+
         mpo[i].insert(shape(0,0,1,2),I);//a

         mpo[i].insert(shape(1,0,1,3),I);//a
         mpo[i].insert(shape(2,1,0,3),I);//a^+

         mpo[i].insert(shape(3,0,0,3),I);
         mpo[i].insert(shape(3,1,1,3),I);

      }

      mpo[L-1].insert(shape(0,0,0,0),O);//zero, you have to put something here: else merge will reduce the rank!

      mpo[L-1].insert(shape(1,0,1,0),I);//a
      mpo[L-1].insert(shape(2,1,0,0),I);//a^+

      mpo[L-1].insert(shape(3,0,0,0),I);
      mpo[L-1].insert(shape(3,1,1,0),I);

      //merge everything together
      TVector<Qshapes<Quantum>,1> qmerge;
      TVector<Dshapes,1> dmerge;

      qmerge[0] = mpo[0].qshape(3);
      dmerge[0] = mpo[0].dshape(3);

      QSTmergeInfo<1> info(qmerge,dmerge);

      QSDArray<4> tmp;
      QSTmerge(mpo[0],info,tmp);

      mpo[0] = tmp;

      for(int i = 1;i < L - 1;++i){

         //first merge the row
         qmerge[0] = mpo[i].qshape(0);
         dmerge[0] = mpo[i].dshape(0);

         info.reset(qmerge,dmerge);

         tmp.clear();

         QSTmerge(info,mpo[i],tmp);

         //then merge the column
         qmerge[0] = tmp.qshape(3);
         dmerge[0] = tmp.dshape(3);

         info.reset(qmerge,dmerge);

         mpo[i].clear();

         QSTmerge(tmp,info,mpo[i]);

      }

      //only merge row for i = L - 1
      qmerge[0] = mpo[L - 1].qshape(0);
      dmerge[0] = mpo[L - 1].dshape(0);

      info.reset(qmerge,dmerge);

      tmp.clear();

      QSTmerge(info,mpo[L - 1],tmp);

      mpo[L - 1] = tmp;

      return mpo;

   }

   /**
    * @return MPO which contains the next to nearest neighbour hopping operator where the nnn-sites are coupled with a factor t_nn
    * i.e. H= -\sum_{i=1}^{L-1} ( a^\dagger_i a_{i+1} +  a^\dagger_{i+1} a_i ) + -t_nn \sum_{i = 1}^{L-2} ( a^\dagger_i a_{i + 2} a^\dagger_{i + 2} a_i )
    */
   MPO nnn_hopping(int L,int d,double t_nn){

      MPO mpo(L);

      //first set the quantumnumbers, before
      Qshapes<Quantum> qp;
      physical(d,qp);

      Qshapes<Quantum> qz; // 0 quantum number
      qz.push_back(Quantum(0));

      Qshapes<Quantum> qi;
      qi.push_back(Quantum::zero());//0
      qi.push_back(Quantum(1));//a
      qi.push_back(Quantum(-1));//a^+
      qi.push_back(Quantum(1));//a
      qi.push_back(Quantum(-1));//a^+
      qi.push_back(Quantum::zero());//I
      
      Qshapes<Quantum> qo;
      qo.push_back(Quantum::zero());//I
      qo.push_back(Quantum(-1));//a^+
      qo.push_back(Quantum(1));//a
      qo.push_back(Quantum(-1));//a^+
      qo.push_back(Quantum(1));//a
      qo.push_back(Quantum::zero());//0

      mpo[0].resize(Quantum::zero(),make_array(qz,qp,-qp,qo));

      //resize & set to 0
      for(int i = 1; i < L-1; ++i)
         mpo[i].resize(Quantum::zero(),make_array(qi,qp,-qp,qo));

      mpo[L-1].resize(Quantum::zero(),make_array(qi,qp,-qp,qz));

      DArray<4> Ip(1,1,1,1);
      Ip = 1;

      DArray<4> Im(1,1,1,1);
      Im = -1;

      DArray<4> Tm(1,1,1,1);
      Tm = -t_nn;

      DArray<4> O(1,1,1,1);
      O = 0;

      //left
      mpo[0].insert(shape(0,0,0,0),Ip);
      mpo[0].insert(shape(0,1,1,0),Ip);

      //hopping
      mpo[0].insert(shape(0,1,0,1),Ip);//a^+
      mpo[0].insert(shape(0,0,1,2),Ip);//a

      //next hopping: nothing on first site
      mpo[0].insert(shape(0,1,0,3),O);//a^+
      mpo[0].insert(shape(0,0,1,4),O);//a

      mpo[0].insert(shape(0,0,0,5),O);//zero, you have to put something here: else merge will reduce the rank!

      //middle
      for(int i = 1;i < L-1;++i){

         mpo[i].insert(shape(0,0,0,0),Ip);
         mpo[i].insert(shape(0,1,1,0),Ip);

         //regular hopping
         mpo[i].insert(shape(0,1,0,1),Ip);//a^+
         mpo[i].insert(shape(0,0,1,2),Ip);//a

         //next-site hopping
         mpo[i].insert(shape(0,1,0,3),Ip);//a^+
         mpo[i].insert(shape(0,0,1,4),Ip);//a

         //fermion sign stuff for the nn-hopping
         mpo[i].insert(shape(1,0,0,3),Ip);
         mpo[i].insert(shape(1,1,1,3),Im);

         mpo[i].insert(shape(2,0,0,4),Ip);
         mpo[i].insert(shape(2,1,1,4),Im);

         //regular hopping input
         mpo[i].insert(shape(1,0,1,5),Im);//a
         mpo[i].insert(shape(2,1,0,5),Im);//a^+

         //next-site hopping input
         mpo[i].insert(shape(3,0,1,5),Tm);//a
         mpo[i].insert(shape(4,1,0,5),Tm);//a^+

         mpo[i].insert(shape(5,0,0,5),Ip);
         mpo[i].insert(shape(5,1,1,5),Ip);

      }

      mpo[L-1].insert(shape(0,0,0,0),O);//zero, you have to put something here: else merge will reduce the rank!

      //regular hopping
      mpo[L-1].insert(shape(1,0,1,0),Im);//a
      mpo[L-1].insert(shape(2,1,0,0),Im);//a^+

      //nn hopping
      mpo[L-1].insert(shape(3,0,1,0),Tm);//a
      mpo[L-1].insert(shape(4,1,0,0),Tm);//a^+

      mpo[L-1].insert(shape(5,0,0,0),Ip);
      mpo[L-1].insert(shape(5,1,1,0),Ip);

      //merge everything together
      TVector<Qshapes<Quantum>,1> qmerge;
      TVector<Dshapes,1> dmerge;

      qmerge[0] = mpo[0].qshape(3);
      dmerge[0] = mpo[0].dshape(3);

      QSTmergeInfo<1> info(qmerge,dmerge);

      QSDArray<4> tmp;
      QSTmerge(mpo[0],info,tmp);

      mpo[0] = tmp;

      for(int i = 1;i < L - 1;++i){

         //first merge the row
         qmerge[0] = mpo[i].qshape(0);
         dmerge[0] = mpo[i].dshape(0);

         info.reset(qmerge,dmerge);

         tmp.clear();

         QSTmerge(info,mpo[i],tmp);

         //then merge the column
         qmerge[0] = tmp.qshape(3);
         dmerge[0] = tmp.dshape(3);

         info.reset(qmerge,dmerge);

         mpo[i].clear();

         QSTmerge(tmp,info,mpo[i]);

      }

      //only merge row for i = L - 1
      qmerge[0] = mpo[L - 1].qshape(0);
      dmerge[0] = mpo[L - 1].dshape(0);

      info.reset(qmerge,dmerge);

      tmp.clear();

      QSTmerge(info,mpo[L - 1],tmp);

      mpo[L - 1] = tmp;

      return mpo;

   }

   /**
    * @return MPO which contains/describes the general one-body operator t_{ij} a^\dagger_i a_j
    * @param t the L x L DArray<2> containing the one-body operator
    */
   MPO one_body(int L,int d,const DArray<2> &t){

      MPO mpo(L);

      //first set the quantumnumbers, before
      Qshapes<Quantum> qp;
      physical(d,qp);

      Qshapes<Quantum> qz; // 0 quantum number
      qz.push_back(Quantum(0));

      Qshapes<Quantum> qi;
      qi.push_back(Quantum::zero());//0

      for(int i = 1;i < L;++i){

         qi.push_back(Quantum(1));//a
         qi.push_back(Quantum(-1));//a^+

      }

      qi.push_back(Quantum::zero());//I


      Qshapes<Quantum> qo;
      qo.push_back(Quantum::zero());//I

      for(int i = 1;i < L;++i){

         qo.push_back(Quantum(-1));//a^+
         qo.push_back(Quantum(1));//a

      }

      qo.push_back(Quantum::zero());//0

      mpo[0].resize(Quantum::zero(),make_array(qz,qp,-qp,qo));

      //resize & set to 0
      for(int i = 1; i < L-1; ++i)
         mpo[i].resize(Quantum::zero(),make_array(qi,qp,-qp,qo));

      mpo[L-1].resize(Quantum::zero(),make_array(qi,qp,-qp,qz));

      DArray<4> Ip(1,1,1,1);
      Ip = 1;

      DArray<4> Im(1,1,1,1);
      Im = -1;

      DArray<4> O(1,1,1,1);
      O = 0;

      //left
      mpo[0].insert(shape(0,0,0,0),Ip);
      mpo[0].insert(shape(0,1,1,0),Ip);

      //nearest neighbour hopping:
      mpo[0].insert(shape(0,1,0,1),Ip);//a^+
      mpo[0].insert(shape(0,0,1,2),Ip);//a

      //long-range hopping: nothing on first site
      for(int i = 1;i < L - 1;++i){

         mpo[0].insert(shape(0,1,0,2*i + 1),O);//a^+
         mpo[0].insert(shape(0,0,1,2*i + 2),O);//a

      }

      DArray<4> T(1,1,1,1);
      T = t(0,0);//diagonal term becomes local

      mpo[0].insert(shape(0,0,0,2*L - 1),T);

      //middle
      for(int i = 1;i < L-1;++i){

         //first row
         mpo[i].insert(shape(0,0,0,0),Ip);
         mpo[i].insert(shape(0,1,1,0),Ip);

         //nearest neighbour hopping:
         mpo[i].insert(shape(0,1,0,1),Ip);//a^+
         mpo[i].insert(shape(0,0,1,2),Ip);//a

         //long-range hopping: nothing on first row
         for(int j = 1;j < L - 1;++j){

            mpo[j].insert(shape(0,1,0,2*j + 1),O);//a^+
            mpo[j].insert(shape(0,0,1,2*j + 2),O);//a

         }

         T = t(i,i);

         mpo[i].insert(shape(0,1,1,2*L - 1),T);//diagonal term becomes local

         //Middle, the fermion signs
         for(int j = 1;j < 2*L -3;++j){

            mpo[i].insert(shape(j,0,0,j + 2),Ip);
            mpo[i].insert(shape(j,1,1,j + 2),Im);

         }

         //last column
         for(int j = 1;j <= i;++j){

            //anni
            T = t(i-j,i);

            mpo[i].insert(shape(2*j - 1,0,1,2*L - 1),T);

            //crea
            T = t(i,i-j);

            mpo[i].insert(shape(2*j,1,0,2*L - 1),T);

         }

         //rest is zero
         for(int j = i + 1;j < L;++j){

            mpo[i].insert(shape(2*j - 1,0,1,2*L - 1),O);

            mpo[i].insert(shape(2*j,1,0,2*L - 1),O);

         }

         //finally the identity in the right bottom
         mpo[i].insert(shape(2*L - 1,0,0,2*L - 1),Ip);
         mpo[i].insert(shape(2*L - 1,1,1,2*L - 1),Ip);

      }

      T = t(L - 1,L - 1);
      mpo[L-1].insert(shape(0,0,0,0),T);

      for(int j = 1;j < L;++j){

         //anni
         T = t(L - 1 - j,L - 1);

         mpo[L - 1].insert(shape(2*j - 1,0,1,0),T);

         //crea
         T = t(L - 1,L - 1 - j);

         mpo[L - 1].insert(shape(2*j,1,0,0),T);

      }

      //finally the identity
      mpo[L - 1].insert(shape(2*L - 1,0,0,0),Ip);
      mpo[L - 1].insert(shape(2*L - 1,1,1,0),Ip);

      //merge everything together
      TVector<Qshapes<Quantum>,1> qmerge;
      TVector<Dshapes,1> dmerge;

      qmerge[0] = mpo[0].qshape(3);
      dmerge[0] = mpo[0].dshape(3);

      QSTmergeInfo<1> info(qmerge,dmerge);

      QSDArray<4> tmp;
      QSTmerge(mpo[0],info,tmp);

      mpo[0] = tmp;

      for(int i = 1;i < L - 1;++i){

         //first merge the row
         qmerge[0] = mpo[i].qshape(0);
         dmerge[0] = mpo[i].dshape(0);

         info.reset(qmerge,dmerge);

         tmp.clear();

         QSTmerge(info,mpo[i],tmp);

         //then merge the column
         qmerge[0] = tmp.qshape(3);
         dmerge[0] = tmp.dshape(3);

         info.reset(qmerge,dmerge);

         mpo[i].clear();

         QSTmerge(tmp,info,mpo[i]);

      }

      //only merge row for i = L - 1
      qmerge[0] = mpo[L - 1].qshape(0);
      dmerge[0] = mpo[L - 1].dshape(0);

      info.reset(qmerge,dmerge);

      tmp.clear();

      QSTmerge(info,mpo[L - 1],tmp);

      mpo[L - 1] = tmp;

      return mpo;

   }

}
