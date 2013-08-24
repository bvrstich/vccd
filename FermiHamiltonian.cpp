#include <iostream>
#include <iomanip>

using std::cout;
using std::endl;
using std::ostream;

#include "include.h"

/**
 * @param qp Qshapes object containing the local quantumnumbers on output, input is destroyed
 */
template<class Q>
void physical(Qshapes<Q> &qp){

   qp.clear();

   qp.push_back(Q(0,0));//zero
   qp.push_back(Q(1,0));//up
   qp.push_back(Q(0,1));//down
   qp.push_back(Q(1,1));//up down pair

}

/**
 * construct an MPO of length L which creates a fermion at site 'site' with spin +1/2 if spin == 0 and spin -1/2 if spin == 1
 */
template<class Q>
MPO<Q> creator(int L,int site,int spin){

   MPO<Q> mpo(L);

   //first set the quantumnumbers, before
   Qshapes<Q> qp;
   physical(qp);

   Qshapes<Q> qz; // 0 quantum number
   qz.push_back(Q(0,0));

   Qshapes<Q> qt; // total quantum number

   if(spin == 0)
      qt.push_back(Q(1,0));
   else
      qt.push_back(Q(0,1));

   //resize & set to 0
   for(int i = 0; i < site; ++i)
      mpo[i].resize(Q::zero(),make_array(qz,qp,-qp,qz));

   //the quantumnumbers on the site of the operation
   mpo[site].resize(Q::zero(),make_array(qz,qp,-qp,-qt));

   //the quantumnumbers after
   for(int i = site + 1;i < L;++i)
      mpo[i].resize(Q::zero(),make_array(qt,qp,-qp,-qt));

   DArray<4> Ip(1,1,1,1);
   Ip = 1;

   DArray<4> Im(1,1,1,1);
   Im = -1;

   //fill it up before: signs
   for(int i = 0;i < site;++i){

      //(-1)^N_site
      mpo[i].insert(shape(0,0,0,0),Ip);
      mpo[i].insert(shape(0,1,1,0),Im);
      mpo[i].insert(shape(0,2,2,0),Im);
      mpo[i].insert(shape(0,3,3,0),Ip);

   }

   //on
   if(spin == 0){//up particle

      mpo[site].insert(shape(0,2,0,0),Ip);
      mpo[site].insert(shape(0,3,1,0),Ip);

   }
   else{//down particle: sign issue

      mpo[site].insert(shape(0,1,0,0),Ip);
      mpo[site].insert(shape(0,3,2,0),Im);

   }

   //after the operator: only identity
   for(int i = site + 1;i < L;++i){

      mpo[i].insert(shape(0,0,0,0),Ip);
      mpo[i].insert(shape(0,1,1,0),Ip);
      mpo[i].insert(shape(0,2,2,0),Ip);
      mpo[i].insert(shape(0,3,3,0),Ip);

   }

   return mpo;

}

/**
 * construct an MPO of length L which annihilates a fermion at site 'site' with spin +1/2 if spin == 0 and spin -1/2 if spin == 1
 */
template<class Q>
MPO<Q> annihilator(int L,int site,int spin){

   MPO<Q> mpo(L);

   //first set the quantumnumbers, before
   Qshapes<Q> qp;
   physical(qp);

   Qshapes<Q> qz; // 0 quantum number
   qz.push_back(Q(0,0));

   Qshapes<Q> qt; // total quantum number

   if(spin == 0)
      qt.push_back(Q(-1,0));
   else
      qt.push_back(Q(0,-1));

   //resize & set to 0
   for(int i = 0; i < site; ++i)
      mpo[i].resize(Q::zero(),make_array(qz,qp,-qp,qz));

   //the quantumnumbers on the site of the operation
   mpo[site].resize(Q::zero(),make_array(qz,qp,-qp,-qt));

   //the quantumnumbers after
   for(int i = site + 1;i < L;++i)
      mpo[i].resize(Q::zero(),make_array(qt,qp,-qp,-qt));

   DArray<4> Ip(1,1,1,1);
   Ip = 1;

   DArray<4> Im(1,1,1,1);
   Im = -1;

   //fill it up before: signs
   for(int i = 0;i < site;++i){

      //(-1)^N_site
      mpo[i].insert(shape(0,0,0,0),Ip);
      mpo[i].insert(shape(0,1,1,0),Im);
      mpo[i].insert(shape(0,2,2,0),Im);
      mpo[i].insert(shape(0,3,3,0),Ip);

   }

   //on
   if(spin == 0){//up particle

      mpo[site].insert(shape(0,0,2,0),Ip);
      mpo[site].insert(shape(0,1,3,0),Ip);

   }
   else{//down particle: sign issue

      mpo[site].insert(shape(0,0,1,0),Ip);
      mpo[site].insert(shape(0,2,3,0),Im);

   }

   //after the operator: only identity
   for(int i = site + 1;i < L;++i){

      mpo[i].insert(shape(0,0,0,0),Ip);
      mpo[i].insert(shape(0,1,1,0),Ip);
      mpo[i].insert(shape(0,2,2,0),Ip);
      mpo[i].insert(shape(0,3,3,0),Ip);

   }

   return mpo;

}

/**
 * construct an MPO which returns the local particle number operator a^+_i a_i
 */
template<class Q>
MPO<Q> n_loc(int L,int site){

   MPO<Q> mpo(L);

   //first set the quantumnumbers, before
   Qshapes<Q> qp;
   physical(qp);

   Qshapes<Q> qz; // 0 quantum number
   qz.push_back(Q::zero());

   //resize & set to 0
   for(int i = 0; i < L; ++i)
      mpo[i].resize(Q::zero(),make_array(qz,qp,-qp,qz));

   DArray<4> I(1,1,1,1);
   I = 1;

   DArray<4> n(1,1,1,1);
   n = 2;

   //fill it up before:unit
   for(int i = 0;i < site;++i){

      mpo[i].insert(shape(0,0,0,0),I);
      mpo[i].insert(shape(0,1,1,0),I);
      mpo[i].insert(shape(0,2,2,0),I);
      mpo[i].insert(shape(0,3,3,0),I);

   }

   //on
   mpo[site].insert(shape(0,1,1,0),I);
   mpo[site].insert(shape(0,2,2,0),I);
   mpo[site].insert(shape(0,3,3,0),n);

   //after the operator
   for(int i = site + 1;i < L;++i){

      mpo[i].insert(shape(0,0,0,0),I);
      mpo[i].insert(shape(0,1,1,0),I);
      mpo[i].insert(shape(0,2,2,0),I);
      mpo[i].insert(shape(0,3,3,0),I);

   }

   return mpo;

}

/**
 * @return MPO which returns the total particle number operator \sum_i a^+_i a_i
 */
template<class Q>
MPO<Q> N_tot(int L){

   MPO<Q> mpo(L);

   //first set the quantumnumbers, before
   Qshapes<Q> qp;
   physical(qp);

   Qshapes<Q> qz; // 0 quantum number
   qz.push_back(Q::zero());

   Qshapes<Q> qi;
   qi.push_back(Q::zero());
   qi.push_back(Q::zero());

   Qshapes<Q> qo;
   qo.push_back(Q::zero());
   qo.push_back(Q::zero());

   mpo[0].resize(Q::zero(),make_array(qz,qp,-qp,qo));

   //resize & set to 0
   for(int i = 1; i < L-1; ++i)
      mpo[i].resize(Q::zero(),make_array(qi,qp,-qp,qo));

   mpo[L-1].resize(Q::zero(),make_array(qi,qp,-qp,qz));

   DArray<4> I(1,1,1,1);
   I = 1;

   DArray<4> n(1,1,1,1);
   n = 2;

   //left
   mpo[0].insert(shape(0,0,0,0),I);
   mpo[0].insert(shape(0,1,1,0),I);
   mpo[0].insert(shape(0,2,2,0),I);
   mpo[0].insert(shape(0,3,3,0),I);

   mpo[0].insert(shape(0,1,1,1),I);
   mpo[0].insert(shape(0,2,2,1),I);
   mpo[0].insert(shape(0,3,3,1),n);

   //middle
   for(int i = 1;i < L-1;++i){

      mpo[i].insert(shape(0,0,0,0),I);
      mpo[i].insert(shape(0,1,1,0),I);
      mpo[i].insert(shape(0,2,2,0),I);
      mpo[i].insert(shape(0,3,3,0),I);

      mpo[i].insert(shape(0,1,1,1),I);
      mpo[i].insert(shape(0,2,2,1),I);
      mpo[i].insert(shape(0,3,3,1),n);

      mpo[i].insert(shape(1,0,0,1),I);
      mpo[i].insert(shape(1,1,1,1),I);
      mpo[i].insert(shape(1,2,2,1),I);
      mpo[i].insert(shape(1,3,3,1),I);

   }

   //right
   mpo[L-1].insert(shape(0,1,1,0),I);
   mpo[L-1].insert(shape(0,2,2,0),I);
   mpo[L-1].insert(shape(0,3,3,0),n);

   mpo[L-1].insert(shape(1,0,0,0),I);
   mpo[L-1].insert(shape(1,1,1,0),I);
   mpo[L-1].insert(shape(1,2,2,0),I);
   mpo[L-1].insert(shape(1,3,3,0),I);

   //merge everything together
   TVector<Qshapes<Q>,1> qmerge;
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
 * @return MPO which returns the total number of spin-up particles
 */
template<class Q>
MPO<Q> n_up_tot(int L){

   MPO<Q> mpo(L);

   //first set the quantumnumbers, before
   Qshapes<Q> qp;
   physical(qp);

   Qshapes<Q> qz; // 0 quantum number
   qz.push_back(Q::zero());

   Qshapes<Q> qi;
   qi.push_back(Q::zero());
   qi.push_back(Q::zero());

   Qshapes<Q> qo;
   qo.push_back(Q::zero());
   qo.push_back(Q::zero());

   mpo[0].resize(Q::zero(),make_array(qz,qp,-qp,qo));

   //resize & set to 0
   for(int i = 1; i < L-1; ++i)
      mpo[i].resize(Q::zero(),make_array(qi,qp,-qp,qo));

   mpo[L-1].resize(Q::zero(),make_array(qi,qp,-qp,qz));

   DArray<4> I(1,1,1,1);
   I = 1;

   //left
   mpo[0].insert(shape(0,0,0,0),I);
   mpo[0].insert(shape(0,1,1,0),I);
   mpo[0].insert(shape(0,2,2,0),I);
   mpo[0].insert(shape(0,3,3,0),I);

   mpo[0].insert(shape(0,2,2,1),I);
   mpo[0].insert(shape(0,3,3,1),I);

   //middle
   for(int i = 1;i < L-1;++i){

      mpo[i].insert(shape(0,0,0,0),I);
      mpo[i].insert(shape(0,1,1,0),I);
      mpo[i].insert(shape(0,2,2,0),I);
      mpo[i].insert(shape(0,3,3,0),I);

      mpo[i].insert(shape(0,2,2,1),I);
      mpo[i].insert(shape(0,3,3,1),I);

      mpo[i].insert(shape(1,0,0,1),I);
      mpo[i].insert(shape(1,1,1,1),I);
      mpo[i].insert(shape(1,2,2,1),I);
      mpo[i].insert(shape(1,3,3,1),I);

   }

   //right
   mpo[L-1].insert(shape(0,2,2,0),I);
   mpo[L-1].insert(shape(0,3,3,0),I);

   mpo[L-1].insert(shape(1,0,0,0),I);
   mpo[L-1].insert(shape(1,1,1,0),I);
   mpo[L-1].insert(shape(1,2,2,0),I);
   mpo[L-1].insert(shape(1,3,3,0),I);

   //merge everything together
   TVector<Qshapes<Q>,1> qmerge;
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
 * @return MPO which returns the total number of spin-down particles
 */
template<class Q>
MPO<Q> n_down_tot(int L){

   MPO<Q> mpo(L);

   //first set the quantumnumbers, before
   Qshapes<Q> qp;
   physical(qp);

   Qshapes<Q> qz; // 0 quantum number
   qz.push_back(Q::zero());

   Qshapes<Q> qi;
   qi.push_back(Q::zero());
   qi.push_back(Q::zero());

   Qshapes<Q> qo;
   qo.push_back(Q::zero());
   qo.push_back(Q::zero());

   mpo[0].resize(Q::zero(),make_array(qz,qp,-qp,qo));

   //resize & set to 0
   for(int i = 1; i < L-1; ++i)
      mpo[i].resize(Q::zero(),make_array(qi,qp,-qp,qo));

   mpo[L-1].resize(Q::zero(),make_array(qi,qp,-qp,qz));

   DArray<4> I(1,1,1,1);
   I = 1;

   //left
   mpo[0].insert(shape(0,0,0,0),I);
   mpo[0].insert(shape(0,1,1,0),I);
   mpo[0].insert(shape(0,2,2,0),I);
   mpo[0].insert(shape(0,3,3,0),I);

   mpo[0].insert(shape(0,1,1,1),I);
   mpo[0].insert(shape(0,3,3,1),I);

   //middle
   for(int i = 1;i < L-1;++i){

      mpo[i].insert(shape(0,0,0,0),I);
      mpo[i].insert(shape(0,1,1,0),I);
      mpo[i].insert(shape(0,2,2,0),I);
      mpo[i].insert(shape(0,3,3,0),I);

      mpo[i].insert(shape(0,1,1,1),I);
      mpo[i].insert(shape(0,3,3,1),I);

      mpo[i].insert(shape(1,0,0,1),I);
      mpo[i].insert(shape(1,1,1,1),I);
      mpo[i].insert(shape(1,2,2,1),I);
      mpo[i].insert(shape(1,3,3,1),I);

   }

   //right
   mpo[L-1].insert(shape(0,1,1,0),I);
   mpo[L-1].insert(shape(0,3,3,0),I);

   mpo[L-1].insert(shape(1,0,0,0),I);
   mpo[L-1].insert(shape(1,1,1,0),I);
   mpo[L-1].insert(shape(1,2,2,0),I);
   mpo[L-1].insert(shape(1,3,3,0),I);

   //merge everything together
   TVector<Qshapes<Q>,1> qmerge;
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
 * @return MPO which contains the hubbard model hamiltonian
 */
template<class Q>
MPO<Q> hubbard(int L,double U){

   MPO<Q> mpo(L);

   //first set the quantumnumbers, before
   Qshapes<Q> qp;
   physical(qp);

   Qshapes<Q> qz; // 0 quantum number
   qz.push_back(Q::zero());

   Qshapes<Q> qi;
   qi.push_back(Q::zero());//on-site rep
   qi.push_back(Q(1,0));//a^+_up
   qi.push_back(Q(-1,0));//a_up
   qi.push_back(Q(0,1));//a^+_down
   qi.push_back(Q(0,-1));//a_up_down
   qi.push_back(Q::zero());//I

   Qshapes<Q> qo;
   qo.push_back(Q::zero());//I
   qo.push_back(Q(-1,0));//a^+_up
   qo.push_back(Q(1,0));//a_up
   qo.push_back(Q(0,-1));//a^+_down
   qo.push_back(Q(0,1));//a_up_down
   qo.push_back(Q::zero());//on-site rep

   mpo[0].resize(Q::zero(),make_array(qz,qp,-qp,qo));

   //resize & set to 0
   for(int i = 1; i < L-1; ++i)
      mpo[i].resize(Q::zero(),make_array(qi,qp,-qp,qo));

   mpo[L-1].resize(Q::zero(),make_array(qi,qp,-qp,qz));

   DArray<4> Ip(1,1,1,1);
   Ip = 1;

   DArray<4> Im(1,1,1,1);
   Im = -1;

   DArray<4> U_op(1,1,1,1);
   U_op = U;

   //left
   //identity
   mpo[0].insert(shape(0,0,0,0),Ip);
   mpo[0].insert(shape(0,1,1,0),Ip);
   mpo[0].insert(shape(0,2,2,0),Ip);
   mpo[0].insert(shape(0,3,3,0),Ip);

   //a^+_up (-1)^{n_down}
   mpo[0].insert(shape(0,2,0,1),Ip);
   mpo[0].insert(shape(0,3,1,1),Im);

   //a_up (-1)^{n_down}
   mpo[0].insert(shape(0,0,2,2),Ip);
   mpo[0].insert(shape(0,1,3,2),Im);

   //a^+_down
   mpo[0].insert(shape(0,1,0,3),Ip);
   mpo[0].insert(shape(0,3,2,3),Ip);

   //a_down
   mpo[0].insert(shape(0,0,1,4),Ip);
   mpo[0].insert(shape(0,2,3,4),Ip);

   mpo[0].insert(shape(0,3,3,5),U_op);//on-site rep

   //middle
   for(int i = 1;i < L-1;++i){

      mpo[i].insert(shape(0,0,0,0),Ip);
      mpo[i].insert(shape(0,1,1,0),Ip);
      mpo[i].insert(shape(0,2,2,0),Ip);
      mpo[i].insert(shape(0,3,3,0),Ip);

      //a^+_up (-1)^{n_down}
      mpo[i].insert(shape(0,2,0,1),Ip);
      mpo[i].insert(shape(0,3,1,1),Im);

      //a_up (-1)^{n_down}
      mpo[i].insert(shape(0,0,2,2),Ip);
      mpo[i].insert(shape(0,1,3,2),Im);

      //a^+_down
      mpo[i].insert(shape(0,1,0,3),Ip);
      mpo[i].insert(shape(0,3,2,3),Ip);

      //a_down
      mpo[i].insert(shape(0,0,1,4),Ip);
      mpo[i].insert(shape(0,2,3,4),Ip);

      mpo[i].insert(shape(0,3,3,5),U_op);//on-site rep

      //a_up
      mpo[i].insert(shape(1,0,2,5),Ip);//on-site rep
      mpo[i].insert(shape(1,1,3,5),Ip);//on-site rep

      //a^+_up
      mpo[i].insert(shape(2,2,0,5),Ip);//on-site rep
      mpo[i].insert(shape(2,3,1,5),Ip);//on-site rep

      //a_down (-1)^n_up
      mpo[i].insert(shape(3,0,1,5),Ip);
      mpo[i].insert(shape(3,2,3,5),Im);

      //a^+_down (-1)^n_up
      mpo[i].insert(shape(4,1,0,5),Ip);
      mpo[i].insert(shape(4,3,2,5),Im);

      //finally Identity
      mpo[i].insert(shape(5,0,0,5),Ip);
      mpo[i].insert(shape(5,1,1,5),Ip);
      mpo[i].insert(shape(5,2,2,5),Ip);
      mpo[i].insert(shape(5,3,3,5),Ip);

   }

   //to the right!
   mpo[L-1].insert(shape(0,3,3,0),U_op);//on-site rep

   //a_up
   mpo[L-1].insert(shape(1,0,2,0),Ip);//on-site rep
   mpo[L-1].insert(shape(1,1,3,0),Ip);//on-site rep

   //a^+_up
   mpo[L-1].insert(shape(2,2,0,0),Ip);//on-site rep
   mpo[L-1].insert(shape(2,3,1,0),Ip);//on-site rep

   //a_down (-1)^n_up
   mpo[L-1].insert(shape(3,0,1,0),Ip);
   mpo[L-1].insert(shape(3,2,3,0),Im);

   //a^+_down (-1)^n_up
   mpo[L-1].insert(shape(4,1,0,0),Ip);
   mpo[L-1].insert(shape(4,3,2,0),Im);

   //finally Identity
   mpo[L-1].insert(shape(5,0,0,0),Ip);
   mpo[L-1].insert(shape(5,1,1,0),Ip);
   mpo[L-1].insert(shape(5,2,2,0),Ip);
   mpo[L-1].insert(shape(5,3,3,0),Ip);

   //merge everything together
   TVector<Qshapes<Q>,1> qmerge;
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
 * @return MPO object of length L containing the T1 operator with coefficients passed through the DArray<2> object t
 */
template<class Q>
MPO<Q> T1(const DArray<2> &t){

   int no = t.shape(0);//number of occupied orbitals
   int nv = t.shape(1);//number of virtual orbitals

   int L = no + nv;

   MPO<Q> mpo(L);

   Qshapes<Q> qp;
   physical(qp);

   Qshapes<Q> qz; // 0 quantum number
   qz.push_back(Q::zero());

   Qshapes<Q> qo;
   qo.push_back(Q::zero());//I
   qo.push_back(Q(1,0));//a_up
   qo.push_back(Q(0,1));//a_down

   DArray<4> Ip(1,1,1,1);
   Ip = 1;

   DArray<4> Im(1,1,1,1);
   Im = -1;

   mpo[0].resize(Q::zero(),make_array(qz,qp,-qp,qo));

   //identity
   mpo[0].insert(shape(0,0,0,0),Ip);
   mpo[0].insert(shape(0,1,1,0),Ip);
   mpo[0].insert(shape(0,2,2,0),Ip);
   mpo[0].insert(shape(0,3,3,0),Ip);

   //a_up (-1)^{n_down}
   mpo[0].insert(shape(0,0,2,1),Ip);
   mpo[0].insert(shape(0,1,3,1),Im);

   //a_down
   mpo[0].insert(shape(0,0,1,2),Ip);
   mpo[0].insert(shape(0,2,3,2),Ip);

   for(int i = 1;i < no - 1;++i){

      qo.clear();

      qo.push_back(Q::zero());//I

      for(int j = 0;j < i + 1;++j){

         qo.push_back(Q(1,0));//a_up
         qo.push_back(Q(0,1));//a_down

      }

      mpo[i].resize(Q::zero(),make_array(-mpo[i - 1].qshape(3),qp,-qp,qo));

      //identity
      mpo[i].insert(shape(0,0,0,0),Ip);
      mpo[i].insert(shape(0,1,1,0),Ip);
      mpo[i].insert(shape(0,2,2,0),Ip);
      mpo[i].insert(shape(0,3,3,0),Ip);

      //a_up (-1)^{n_down}
      mpo[i].insert(shape(0,0,2,1),Ip);
      mpo[i].insert(shape(0,1,3,1),Im);

      //a_down
      mpo[i].insert(shape(0,0,1,2),Ip);
      mpo[i].insert(shape(0,2,3,2),Ip);

      //signs!
      for(int j = 1;j < 2*i + 1;++j){

         //fermion sign!
         mpo[i].insert(shape(j,0,0,j + 2),Ip);
         mpo[i].insert(shape(j,1,1,j + 2),Im);
         mpo[i].insert(shape(j,2,2,j + 2),Im);
         mpo[i].insert(shape(j,3,3,j + 2),Ip);

      }

   }

   //last occupied
   qo.clear();

   for(int j = 0;j < no;++j){

      qo.push_back(Q(1,0));//a_up
      qo.push_back(Q(0,1));//a_down

   }

   mpo[no - 1].resize(Q::zero(),make_array(-mpo[no - 2].qshape(3),qp,-qp,qo));

   //a_up (-1)^{n_down}
   mpo[no - 1].insert(shape(0,0,2,0),Ip);
   mpo[no - 1].insert(shape(0,1,3,0),Im);

   //a_down
   mpo[no - 1].insert(shape(0,0,1,1),Ip);
   mpo[no - 1].insert(shape(0,2,3,1),Ip);

   //signs!
   for(int j = 1;j < 2*no - 1;++j){

      //fermion sign!
      mpo[no - 1].insert(shape(j,0,0,j + 1),Ip);
      mpo[no - 1].insert(shape(j,1,1,j + 1),Im);
      mpo[no - 1].insert(shape(j,2,2,j + 1),Im);
      mpo[no - 1].insert(shape(j,3,3,j + 1),Ip);

   }

   //first virtual
   qo.clear();

   for(int j = 0;j < no;++j){

      qo.push_back(Q(1,0));//a_up
      qo.push_back(Q(0,1));//a_down

   }

   qo.push_back(Q::zero());//last column

   mpo[no].resize(Q::zero(),make_array(-mpo[no - 1].qshape(3),qp,-qp,qo));

   //signs!
   for(int j = 0;j < 2*no;++j){

      //fermion sign!
      mpo[no].insert(shape(j,0,0,j),Ip);
      mpo[no].insert(shape(j,1,1,j),Im);
      mpo[no].insert(shape(j,2,2,j),Im);
      mpo[no].insert(shape(j,3,3,j),Ip);

   }

   DArray<4> Tp(1,1,1,1);
   DArray<4> Tm(1,1,1,1);

   //last column: create!
   for(int j = 0;j < no;++j){

      Tp = t(no - 1 - j,0);
      Tm = -t(no - 1 - j,0);

      //a^+ up
      mpo[no].insert(shape(2*j,2,0,2*no),Tp);
      mpo[no].insert(shape(2*j,3,1,2*no),Tp);

      //a^+ down (-1)^n_up
      mpo[no].insert(shape(2*j + 1,1,0,2*no),Tp);
      mpo[no].insert(shape(2*j + 1,3,2,2*no),Tm);

   }

   //other virtuals
   for(int i = no + 1;i < L - 1;++i){

      mpo[i].resize(Q::zero(),make_array(-qo,qp,-qp,qo));

      //signs!
      for(int j = 0;j < 2*no;++j){

         //fermion sign!
         mpo[i].insert(shape(j,0,0,j),Ip);
         mpo[i].insert(shape(j,1,1,j),Im);
         mpo[i].insert(shape(j,2,2,j),Im);
         mpo[i].insert(shape(j,3,3,j),Ip);

      }

      //last column: create!
      for(int j = 0;j < no;++j){

         Tp = t(no - 1 - j,i - no);
         Tm = -t(no - 1 - j,i - no);

         //a^+ up
         mpo[i].insert(shape(2*j,2,0,2*no),Tp);
         mpo[i].insert(shape(2*j,3,1,2*no),Tp);

         //a^+ down (-1)^n_up
         mpo[i].insert(shape(2*j + 1,1,0,2*no),Tp);
         mpo[i].insert(shape(2*j + 1,3,2,2*no),Tm);

      }

      //last element:identity
      mpo[i].insert(shape(2*no,0,0,2*no),Ip);
      mpo[i].insert(shape(2*no,1,1,2*no),Ip);
      mpo[i].insert(shape(2*no,2,2,2*no),Ip);
      mpo[i].insert(shape(2*no,3,3,2*no),Ip);

   }

   //finally the last virtual
   mpo[L-1].resize(Q::zero(),make_array(-qo,qp,-qp,qz));

   for(int j = 0;j < no;++j){

      Tp = t(no - 1 - j,nv - 1);
      Tm = -t(no - 1 - j,nv - 1);

      //a^+ up
      mpo[L-1].insert(shape(2*j,2,0,0),Tp);
      mpo[L-1].insert(shape(2*j,3,1,0),Tp);

      //a^+ down (-1)^n_up
      mpo[L-1].insert(shape(2*j + 1,1,0,0),Tp);
      mpo[L-1].insert(shape(2*j + 1,3,2,0),Tm);

   }

   //last element:identity
   mpo[L-1].insert(shape(2*no,0,0,0),Ip);
   mpo[L-1].insert(shape(2*no,1,1,0),Ip);
   mpo[L-1].insert(shape(2*no,2,2,0),Ip);
   mpo[L-1].insert(shape(2*no,3,3,0),Ip);

   //merge everything together
   TVector<Qshapes<Q>,1> qmerge;
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
 * @return MPO object of length L containing the T2 operator with coefficients passed through the DArray<4> object t
 */
template<class Q>
MPO<Q> T2(const DArray<4> &t){

   int no = t.shape(0);//number of occupied orbitals
   int nv = t.shape(2);//number of virtual orbitals

   int L = no + nv;

   MPO<Q> mpo(L);

   Qshapes<Q> qp;
   physical(qp);

   Qshapes<Q> qz; // 0 quantum number
   qz.push_back(Q::zero());

   Qshapes<Q> qo;
   qo.push_back(Q::zero());//I
   qo.push_back(Q(0,1));//a_down
   qo.push_back(Q(1,0));//a_up
   qo.push_back(Q(1,1));//a_down a_up

   DArray<4> Ip(1,1,1,1);
   Ip = 1;

   DArray<4> Im(1,1,1,1);
   Im = -1;

   mpo[0].resize(Q::zero(),make_array(qz,qp,-qp,qo));

   std::vector<int> v(1);
   v[0] = 1;

   std::vector< std::vector<int> > ostates;

   //identity
   mpo[0].insert(shape(0,0,0,0),Ip);
   mpo[0].insert(shape(0,1,1,0),Ip);
   mpo[0].insert(shape(0,2,2,0),Ip);
   mpo[0].insert(shape(0,3,3,0),Ip);

   ostates.push_back(v);

   //a_down
   mpo[0].insert(shape(0,0,2,1),Ip);
   mpo[0].insert(shape(0,1,3,1),Ip);

   v.resize(2);
   v[0] = 0;//site
   v[1] = 1;//spin down

   ostates.push_back(v);

   //a_up (-1)^n_down
   mpo[0].insert(shape(0,0,1,2),Ip);
   mpo[0].insert(shape(0,2,3,2),Im);

   v[0] = 0;//site
   v[1] = 0;//spin up

   ostates.push_back(v);

   //a_down a_up
   mpo[0].insert(shape(0,0,3,3),Ip);

   v.resize(4);

   v[0] = 0;//site
   v[1] = 1;//down
   v[2] = 0;//site
   v[3] = 0;//up

   ostates.push_back(v);

   std::vector< std::vector<int> > istates;
   istates = ostates;

   for(int i = 1;i < no - 1;++i){

      ostates.clear();

      qo.clear();

      qo.push_back(Q::zero());//I

      ostates.push_back(istates[0]);

      qo.push_back(Q(0,1));//a_down
      qo.push_back(Q(1,0));//a_up

      v.resize(2);
      v[0] = i;//site i
      v[1] = 1;//spin down

      ostates.push_back(v);

      v[0] = i;//site i
      v[1] = 0;//spin up

      ostates.push_back(v);

      for(int j = 1;j < i + 1;++j){

         qo.push_back(Q(0,1));//a_down
         ostates.push_back(istates[2*j - 1]);

         qo.push_back(Q(1,0));//a_up
         ostates.push_back(istates[2*j]);

      }

      qo.push_back(Q(1,1));//a_up a_down

      v.resize(4);
      v[0] = i;
      v[1] = 1;
      v[2] = i;
      v[3] = 0;

      ostates.push_back(v);

      //add a down
      for(int j = 1;j < i + 1;++j){

         qo.push_back(Q(0,2));//a_down a_down

         v[0] = i;
         v[1] = 1;
         v[2] = istates[2*j - 1][0];
         v[3] = istates[2*j - 1][1];

         ostates.push_back(v);

         qo.push_back(Q(1,1));//a_up a_down

         v[0] = i;
         v[1] = 1;
         v[2] = istates[2*j][0];
         v[3] = istates[2*j][1];

         ostates.push_back(v);

      }

      //add an up 
      for(int j = 1;j < i + 1;++j){

         qo.push_back(Q(1,1));//a_down a_up

         v[0] = i;
         v[1] = 0;
         v[2] = istates[2*j - 1][0];
         v[3] = istates[2*j - 1][1];

         ostates.push_back(v);

         qo.push_back(Q(2,0));//a_up a_up

         v[0] = i;
         v[1] = 0;
         v[2] = istates[2*j][0];
         v[3] = istates[2*j][1];

         ostates.push_back(v);

      }

      //identity for the pairs
      for(int j = i - 1;j > 0;--j){

         //two on site
         qo.push_back(Q(1,1));//a_up a_down

         for(int k = j - 1;k >= 0;--k){

            qo.push_back(Q(0,2));//a_down a_down
            qo.push_back(Q(1,1));//a_up a_down

         }

         for(int k = j - 1;k >= 0;--k){

            qo.push_back(Q(1,1));//a_up a_down
            qo.push_back(Q(2,0));//a_up a_up

         }

      }

      //a_0_up a_0_down
      qo.push_back(Q(1,1));

      for(int j = 2*i + 1;j < istates.size();++j)
         ostates.push_back(istates[j]);

      mpo[i].resize(Q::zero(),make_array(-mpo[i - 1].qshape(3),qp,-qp,qo));

      int column = 0;

      //identity
      mpo[i].insert(shape(0,0,0,0),Ip);
      mpo[i].insert(shape(0,1,1,0),Ip);
      mpo[i].insert(shape(0,2,2,0),Ip);
      mpo[i].insert(shape(0,3,3,0),Ip);

      //a_down
      mpo[i].insert(shape(0,0,2,1),Ip);
      mpo[i].insert(shape(0,1,3,1),Ip);

      //a_up (-1)^{n_down}
      mpo[i].insert(shape(0,0,1,2),Ip);
      mpo[i].insert(shape(0,2,3,2),Im);

      column = 3;

      //signs!
      for(int j = 1;j < 2*i + 1;++j){

         //fermion sign!
         mpo[i].insert(shape(j,0,0,column),Ip);
         mpo[i].insert(shape(j,1,1,column),Im);
         mpo[i].insert(shape(j,2,2,column),Im);
         mpo[i].insert(shape(j,3,3,column),Ip);

         column++;

      }

      //double remove
      mpo[i].insert(shape(0,0,3,column),Ip);
      column++;

      //single remove down (-1)^n_up
      for(int j = 1;j < 2*i + 1;++j){

         //fermion sign!
         mpo[i].insert(shape(j,0,2,column),Ip);
         mpo[i].insert(shape(j,1,3,column),Im);
         column++;

      }

      //single remove up
      for(int j = 1;j < 2*i + 1;++j){

         //fermion sign!
         mpo[i].insert(shape(j,0,1,column),Ip);
         mpo[i].insert(shape(j,2,3,column),Ip);
         ++column;

      }

      int row = 2*i + 1;

      for(int j = 0;j < i*(2*i - 1);++j){

         //finally identity
         mpo[i].insert(shape(row,0,0,column),Ip);
         mpo[i].insert(shape(row,1,1,column),Ip);
         mpo[i].insert(shape(row,2,2,column),Ip);
         mpo[i].insert(shape(row,3,3,column),Ip);
         ++row;
         ++column;

      }

      istates = ostates;

   }

   ostates.clear();

   //last occupied: i = no - 1
   qo.clear();

   qo.push_back(Q(1,1));//a_up a_down

   v.resize(4);
   v[0] = no - 1;
   v[1] = 1;
   v[2] = no - 1;
   v[3] = 0;

   ostates.push_back(v);

   for(int j = 1;j < no;++j){

      qo.push_back(Q(0,2));//a_down a_down

      v[0] = no - 1;
      v[1] = 1;
      v[2] = istates[2*j - 1][0];
      v[3] = istates[2*j - 1][1];

      ostates.push_back(v);

      qo.push_back(Q(1,1));//a_up a_down

      v[0] = no - 1;
      v[1] = 1;
      v[2] = istates[2*j][0];
      v[3] = istates[2*j][1];

      ostates.push_back(v);

   }

   for(int j = 1;j < no;++j){

      qo.push_back(Q(1,1));//a_down a_up

      v[0] = no - 1;
      v[1] = 0;
      v[2] = istates[2*j - 1][0];
      v[3] = istates[2*j - 1][1];

      ostates.push_back(v);

      qo.push_back(Q(2,0));//a_up a_up

      v[0] = no - 1;
      v[1] = 0;
      v[2] = istates[2*j][0];
      v[3] = istates[2*j][1];

      ostates.push_back(v);

   }

   //identity for the pairs from previous iteration
   for(int j = no - 2;j > 0;--j){

      //two on site
      qo.push_back(Q(1,1));//a_up a_down

      for(int k = j - 1;k >= 0;--k){

         qo.push_back(Q(0,2));//a_down a_down
         qo.push_back(Q(1,1));//a_down a_up

      }

      for(int k = j - 1;k >= 0;--k){

         qo.push_back(Q(1,1));//a_up a_down
         qo.push_back(Q(2,0));//a_up a_up

      } 

   }

   //a_0_up a_0_down
   qo.push_back(Q(1,1));

   for(int j = 2*no - 1;j < istates.size();++j)
      ostates.push_back(istates[j]);

   //start filling
   mpo[no - 1].resize(Q::zero(),make_array(-mpo[no - 2].qshape(3),qp,-qp,qo));

   int column = 0;

   //double remove
   mpo[no - 1].insert(shape(0,0,3,column),Ip);
   column++;

   //single remove down
   for(int j = 1;j < 2*no - 1;++j){

      //fermion sign!
      mpo[no - 1].insert(shape(j,0,2,column),Ip);
      mpo[no - 1].insert(shape(j,1,3,column),Im);
      column++;

   }

   //single remove up
   for(int j = 1;j < 2*no - 1;++j){

      mpo[no - 1].insert(shape(j,0,1,column),Ip);
      mpo[no - 1].insert(shape(j,2,3,column),Ip);
      ++column;

   }

   int row = 2*no - 1;

   for(int j = 0;j < (no - 1)*(2*(no - 1) - 1);++j){

      //finally identity
      mpo[no - 1].insert(shape(row,0,0,column),Ip);
      mpo[no - 1].insert(shape(row,1,1,column),Ip);
      mpo[no - 1].insert(shape(row,2,2,column),Ip);
      mpo[no - 1].insert(shape(row,3,3,column),Ip);

      ++row;
      ++column;

   }

   istates = ostates;

   ostates.clear();

   //first virtual
   qo.clear();

   qo.push_back(Q::zero());//first column completely closed 

   v.resize(1);
   v[0] = 1;
   ostates.push_back(v);

   //first all the virtuals
   for(int j = no + 1;j < L;++j){

      qo.push_back(Q(0,1));//go to 1 down removed

      v.resize(2);
      v[0] = j;
      v[1] = 1;
      ostates.push_back(v);

      qo.push_back(Q(1,0));//go to 1 up removed

      v[0] = j;
      v[1] = 0;
      ostates.push_back(v);

   }

   //the rest is identity
   for(int i = 0;i < istates.size();++i){

      ostates.push_back(istates[i]);
      qo.push_back(mpo[no - 1].qshape(3)[i]);

   }

   //start the filling
   mpo[no].resize(Q::zero(),make_array(-mpo[no - 1].qshape(3),qp,-qp,qo));

   //first the local term
   for(int j = 0;j < istates.size();++j){

      if(mpo[no].qshape(0)[j].gn_up() == -1 && mpo[no].qshape(0)[j].gn_down() == -1){//insert double create

         DArray<4> Tp(1,1,1,1);

         if(istates[j][1] == 1)//down spin is first
            Tp = t(istates[j][2],istates[j][0],0,0);
         else//up spin is first
            Tp = -t(istates[j][2],istates[j][0],0,0);

         //double create: first column
         mpo[no].insert(shape(j,3,0,0),Tp);

      }

   }

   //then the rest of the virtuals
   for(int j = no + 1;j < L;++j){

      int vjnd = j - no;//virtual index

      //first to one down removed
      for(int k = 0;k < istates.size();++k){//these are the rows, sum over pairs

         if(mpo[no].qshape(0)[k].gn_up() == -1 && mpo[no].qshape(0)[k].gn_down() == -1){//insert create up, so from (-1,-1) -> (0,1)

            DArray<4> Tp(1,1,1,1);
            DArray<4> Tm(1,1,1,1);

            if(istates[k][2] == istates[k][0]){//is it a double occupied site coming in?

               Tp = t(istates[k][2],istates[k][0],0,vjnd);
               Tm = -t(istates[k][2],istates[k][0],0,vjnd);

            }
            else{//two different sites coming in: up down or down up

               if(istates[k][1] == 1){// down up

                  Tp = t(istates[k][2],istates[k][0],0,vjnd);
                  Tm = -t(istates[k][2],istates[k][0],0,vjnd);

               }
               else{//up down

                  Tp = -t(istates[k][2],istates[k][0],vjnd,0);
                  Tm = t(istates[k][2],istates[k][0],vjnd,0);

               }

            }

            //create up: 
            mpo[no].insert(shape(k,1,0,2*vjnd-1),Tp);
            mpo[no].insert(shape(k,3,2,2*vjnd-1),Tm);

         }
         else if(mpo[no].qshape(0)[k].gn_up() == 0 && mpo[no].qshape(0)[k].gn_down() == -2){//insert create down: from (0,2) -> (0,-1)

            DArray<4> Tp(1,1,1,1);

            Tp = t(istates[k][2],istates[k][0],0,vjnd) - t(istates[k][2],istates[k][0],vjnd,0);

            //create down: 
            mpo[no].insert(shape(k,2,0,2*vjnd-1),Tp);
            mpo[no].insert(shape(k,3,1,2*vjnd-1),Tp);

         }

      }

      //then to one up removed
      for(int k = 0;k < istates.size();++k){

         if(mpo[no].qshape(0)[k].gn_up() == -1 && mpo[no].qshape(0)[k].gn_down() == -1){//insert create down, from (1,1) -> (-1,0)

            DArray<4> Tp(1,1,1,1);

            if(istates[k][2] == istates[k][0])//is it a double occupied site coming in?
               Tp = -t(istates[k][2],istates[k][0],0,vjnd);
            else{//two different sites coming in: up down or down up

               if(istates[k][1] == 1)// down up
                  Tp = -t(istates[k][2],istates[k][0],vjnd,0);
               else//up down
                  Tp = t(istates[k][2],istates[k][0],0,vjnd);

            }

            //create down: 
            mpo[no].insert(shape(k,2,0,2*vjnd),Tp);
            mpo[no].insert(shape(k,3,1,2*vjnd),Tp);

         }
         else if(mpo[no].qshape(0)[k].gn_up() == -2 && mpo[no].qshape(0)[k].gn_down() == 0){//insert create up : from (2,0) -> (-1,0)

            DArray<4> Tp(1,1,1,1);
            DArray<4> Tm(1,1,1,1);

            Tp = t(istates[k][2],istates[k][0],0,vjnd) - t(istates[k][2],istates[k][0],vjnd,0);
            Tm = t(istates[k][2],istates[k][0],vjnd,0) + t(istates[k][2],istates[k][0],0,vjnd);

            //create up:
            mpo[no].insert(shape(k,1,0,2*vjnd),Tp);
            mpo[no].insert(shape(k,3,2,2*vjnd),Tm);

         }

      }

   }

   //fill the remaining part unity
   for(int j = 0;j < istates.size();++j){

      mpo[no].insert(shape(j,0,0,j + 2*nv - 1),Ip);
      mpo[no].insert(shape(j,1,1,j + 2*nv - 1),Ip);
      mpo[no].insert(shape(j,2,2,j + 2*nv - 1),Ip);
      mpo[no].insert(shape(j,3,3,j + 2*nv - 1),Ip);

   }

   //this keeps track of all the pairs
   std::vector< std::vector<int> > pairs = istates;
   istates = ostates;

   //next all but the last virtuals:
   for(int i = no + 1;i < L - 1;++i){

      //current virtual index
      int vind = i - no;

      ostates.clear();

      qo.clear();

      qo.push_back(Q::zero());//first column completely closed 

      v.resize(1);
      v[0] = 1;
      ostates.push_back(v);

      //first all the remaining the virtuals
      for(int j = i + 1;j < L;++j){

         qo.push_back(Q(0,1));//go to 1 down removed

         v.resize(2);
         v[0] = j;
         v[1] = 1;
         ostates.push_back(v);

         qo.push_back(Q(1,0));//go to 1 down removed

         v[0] = j;
         v[1] = 0;
         ostates.push_back(v);

      }

      //then just the pairs
      for(int j = 0;j < pairs.size();++j){

         ostates.push_back(pairs[j]);
         qo.push_back(mpo[i - 1].qshape(3)[2*(nv - vind) + 1 + j]);

      }

      //start filling
      mpo[i].resize(Q::zero(),make_array(-mpo[i - 1].qshape(3),qp,-qp,qo));

      //first column is closed
      mpo[i].insert(shape(0,0,0,0),Ip);
      mpo[i].insert(shape(0,1,1,0),Ip);
      mpo[i].insert(shape(0,2,2,0),Ip);
      mpo[i].insert(shape(0,3,3,0),Ip);

      //create a down particle
      mpo[i].insert(shape(1,2,0,0),Ip);
      mpo[i].insert(shape(1,3,1,0),Im);

      //create an up particle
      mpo[i].insert(shape(2,1,0,0),Ip);
      mpo[i].insert(shape(2,3,2,0),Ip);

      row = 2*(nv - vind) + 1;

      //insert the closing pairs
      for(int j = 0;j < pairs.size();++j){

         if(mpo[i].qshape(0)[row].gn_up() == -1 && mpo[i].qshape(0)[row].gn_down() == -1){//insert double create

            DArray<4> Tp(1,1,1,1);

            if(pairs[j][1] == 1)//down spin is first
               Tp = t(pairs[j][2],pairs[j][0],vind,vind);
            else//up spin is first
               Tp = -t(pairs[j][2],pairs[j][0],vind,vind);

            //double create: first column
            mpo[i].insert(shape(row,3,0,0),Tp);

         }

         ++row;

      }

      //insert fermion signs
      for(int j = 3;j < 2*(nv - vind) + 1;++j){

         mpo[i].insert(shape(j,0,0,j-2),Ip);
         mpo[i].insert(shape(j,1,1,j-2),Im);
         mpo[i].insert(shape(j,2,2,j-2),Im);
         mpo[i].insert(shape(j,3,3,j-2),Ip);

      }

      //insert the t coefficients
      for(int j = i + 1;j < L;++j){

         int vjnd = j - no;
         int col = j - i;

         row = 2*(nv - vind) + 1;

         //first to one down removed
         for(int k = 0;k < pairs.size();++k){//these are the rows, sum over pairs

            if(mpo[i].qshape(0)[row].gn_up() == -1 && mpo[i].qshape(0)[row].gn_down() == -1){//insert create up, so from (-1,-1) -> (0,1)

               DArray<4> Tp(1,1,1,1);
               DArray<4> Tm(1,1,1,1);

               if(pairs[k][2] == pairs[k][0]){//is it a double occupied site coming in?

                  Tp = t(pairs[k][2],pairs[k][0],vind,vjnd);
                  Tm = -t(pairs[k][2],pairs[k][0],vind,vjnd);

               }
               else{//two different sites coming in: up down or down up

                  if(pairs[k][1] == 1){// down up

                     Tp = t(pairs[k][2],pairs[k][0],vind,vjnd);
                     Tm = -t(pairs[k][2],pairs[k][0],vind,vjnd);

                  }
                  else{//up down

                     Tp = -t(pairs[k][2],pairs[k][0],vjnd,vind);
                     Tm = t(pairs[k][2],pairs[k][0],vjnd,vind);

                  }

               }

               //create up: 
               mpo[i].insert(shape(row,1,0,2*col-1),Tp);
               mpo[i].insert(shape(row,3,2,2*col-1),Tm);

            }
            else if(mpo[i].qshape(0)[row].gn_up() == 0 && mpo[i].qshape(0)[row].gn_down() == -2){//insert create down: from (0,2) -> (0,-1)

               DArray<4> Tp(1,1,1,1);

               Tp = t(pairs[k][2],pairs[k][0],vind,vjnd) - t(pairs[k][2],pairs[k][0],vjnd,vind);

               //create down: 
               mpo[i].insert(shape(row,2,0,2*col-1),Tp);
               mpo[i].insert(shape(row,3,1,2*col-1),Tp);

            }

            ++row;

         }

         row = 2*(nv - vind) + 1;

         //then to one up removed
         for(int k = 0;k < pairs.size();++k){

            if(mpo[i].qshape(0)[row].gn_up() == -1 && mpo[i].qshape(0)[row].gn_down() == -1){//insert create down, from (1,1) -> (-1,0)

               DArray<4> Tp(1,1,1,1);

               if(pairs[k][2] == pairs[k][0])//is it a double occupied site coming in?
                  Tp = -t(pairs[k][2],pairs[k][0],vind,vjnd);
               else{//two different sites coming in: up down or down up

                  if(pairs[k][1] == 1)// down up
                     Tp = -t(pairs[k][2],pairs[k][0],vjnd,vind);
                  else//up down
                     Tp = t(pairs[k][2],pairs[k][0],vind,vjnd);

               }

               //create down: 
               mpo[i].insert(shape(row,2,0,2*col),Tp);
               mpo[i].insert(shape(row,3,1,2*col),Tp);

            }
            else if(mpo[i].qshape(0)[row].gn_up() == -2 && mpo[i].qshape(0)[row].gn_down() == 0){//insert create up : from (2,0) -> (-1,0)

               DArray<4> Tp(1,1,1,1);
               DArray<4> Tm(1,1,1,1);

               Tp = t(pairs[k][2],pairs[k][0],vind,vjnd) - t(pairs[k][2],pairs[k][0],vjnd,vind);
               Tm = t(pairs[k][2],pairs[k][0],vjnd,vind) + t(pairs[k][2],pairs[k][0],vind,vjnd);

               //create up:
               mpo[i].insert(shape(row,1,0,2*col),Tp);
               mpo[i].insert(shape(row,3,2,2*col),Tm);

            }

            ++row;

         }

      }

      //fill the lower right block with unity
      for(int j = 2*(nv - vind) + 1;j < istates.size();++j){

         mpo[i].insert(shape(j,0,0,j - 2),Ip);
         mpo[i].insert(shape(j,1,1,j - 2),Ip);
         mpo[i].insert(shape(j,2,2,j - 2),Ip);
         mpo[i].insert(shape(j,3,3,j - 2),Ip);

      }

      istates = ostates;

   }

   //last virtual L - 1, closing everything down
   ostates.clear();
   qo.clear();

   qo.push_back(Q::zero());//first column completely closed 

   v.resize(1);
   v[0] = 1;
   ostates.push_back(v);

   mpo[L - 1].resize(Q::zero(),make_array(-mpo[L - 2].qshape(3),qp,-qp,qo));

   //first column is closed
   mpo[L - 1].insert(shape(0,0,0,0),Ip);
   mpo[L - 1].insert(shape(0,1,1,0),Ip);
   mpo[L - 1].insert(shape(0,2,2,0),Ip);
   mpo[L - 1].insert(shape(0,3,3,0),Ip);

   //create a down particle
   mpo[L - 1].insert(shape(1,2,0,0),Ip);
   mpo[L - 1].insert(shape(1,3,1,0),Im);

   //create an up particle
   mpo[L - 1].insert(shape(2,1,0,0),Ip);
   mpo[L - 1].insert(shape(2,3,2,0),Ip);

   row = 3;

   //insert the closing pairs
   for(int j = 0;j < pairs.size();++j){

      if(mpo[L - 1].qshape(0)[row].gn_up() == -1 && mpo[L - 1].qshape(0)[row].gn_down() == -1){//insert double create

         DArray<4> Tp(1,1,1,1);

         if(pairs[j][1] == 1)//down spin is first
            Tp = t(pairs[j][2],pairs[j][0],nv - 1,nv - 1);
         else//up spin is first
            Tp = -t(pairs[j][2],pairs[j][0],nv - 1,nv - 1);

         //double create: first column
         mpo[L - 1].insert(shape(row,3,0,0),Tp);

      }

      ++row;

   }

   //merge everything together
   TVector<Qshapes<Q>,1> qmerge;
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
 * @return MPO object of length L containing a general one-body operator \sum_ab t_ab a^+_as a_bs
 */
template<class Q>
MPO<Q> one_body(const DArray<2> &t){

   int L = t.shape(0);//number of occupied orbitals

   MPO<Q> mpo(L);

   Qshapes<Q> qp;
   physical(qp);

   Qshapes<Q> qz; // 0 quantum number
   qz.push_back(Q::zero());

   Qshapes<Q> qo;
   qo.push_back(Q::zero());//I

   qo.push_back(Q(-1,0));//a^+_up
   qo.push_back(Q(0,-1));//a^+_down

   qo.push_back(Q(1,0));//a_up
   qo.push_back(Q(0,1));//a_down

   qo.push_back(Q(0,0));//a^+_up a_up + a^+_down a_down

   DArray<4> Ip(1,1,1,1);
   Ip = 1;

   DArray<4> Im(1,1,1,1);
   Im = -1;

   mpo[0].resize(Q::zero(),make_array(qz,qp,-qp,qo));

   std::vector<int> v(1);
   v[0] = 1;

   std::vector< std::vector<int> > ostates;

   //identity
   mpo[0].insert(shape(0,0,0,0),Ip);
   mpo[0].insert(shape(0,1,1,0),Ip);
   mpo[0].insert(shape(0,2,2,0),Ip);
   mpo[0].insert(shape(0,3,3,0),Ip);

   ostates.push_back(v);

   //a^+_up (-1)^n_down
   mpo[0].insert(shape(0,1,0,1),Ip);
   mpo[0].insert(shape(0,3,2,1),Im);

   v.resize(3);
   v[0] = 0;//site
   v[1] = 0;//spin up
   v[2] = 0;//creator

   ostates.push_back(v);

   //a^dagger_down 
   mpo[0].insert(shape(0,2,0,2),Ip);
   mpo[0].insert(shape(0,3,1,2),Ip);

   v[0] = 0;//site
   v[1] = 1;//spin down
   v[2] = 0;//creator

   ostates.push_back(v);

   //a_up (-1)^n_down
   mpo[0].insert(shape(0,0,1,2),Ip);
   mpo[0].insert(shape(0,2,3,2),Im);

   v[0] = 0;//site
   v[1] = 0;//spin up
   v[2] = 1;//annihilator

   ostates.push_back(v);

   //a_down
   mpo[0].insert(shape(0,0,2,2),Ip);
   mpo[0].insert(shape(0,1,3,2),Ip);

   v[0] = 0;//site
   v[1] = 1;//spin down
   v[2] = 1;//annihilator

   ostates.push_back(v);
/*
   std::vector< std::vector<int> > istates;
   istates = ostates;

   for(int i = 1;i < no - 1;++i){

      ostates.clear();

      qo.clear();

      qo.push_back(Q::zero());//I

      ostates.push_back(istates[0]);

      qo.push_back(Q(0,1));//a_down
      qo.push_back(Q(1,0));//a_up

      v.resize(2);
      v[0] = i;//site i
      v[1] = 1;//spin down

      ostates.push_back(v);

      v[0] = i;//site i
      v[1] = 0;//spin up

      ostates.push_back(v);

      for(int j = 1;j < i + 1;++j){

         qo.push_back(Q(0,1));//a_down
         ostates.push_back(istates[2*j - 1]);

         qo.push_back(Q(1,0));//a_up
         ostates.push_back(istates[2*j]);

      }

      qo.push_back(Q(1,1));//a_up a_down

      v.resize(4);
      v[0] = i;
      v[1] = 1;
      v[2] = i;
      v[3] = 0;

      ostates.push_back(v);

      //add a down
      for(int j = 1;j < i + 1;++j){

         qo.push_back(Q(0,2));//a_down a_down

         v[0] = i;
         v[1] = 1;
         v[2] = istates[2*j - 1][0];
         v[3] = istates[2*j - 1][1];

         ostates.push_back(v);

         qo.push_back(Q(1,1));//a_up a_down

         v[0] = i;
         v[1] = 1;
         v[2] = istates[2*j][0];
         v[3] = istates[2*j][1];

         ostates.push_back(v);

      }

      //add an up 
      for(int j = 1;j < i + 1;++j){

         qo.push_back(Q(1,1));//a_down a_up

         v[0] = i;
         v[1] = 0;
         v[2] = istates[2*j - 1][0];
         v[3] = istates[2*j - 1][1];

         ostates.push_back(v);

         qo.push_back(Q(2,0));//a_up a_up

         v[0] = i;
         v[1] = 0;
         v[2] = istates[2*j][0];
         v[3] = istates[2*j][1];

         ostates.push_back(v);

      }

      //identity for the pairs
      for(int j = i - 1;j > 0;--j){

         //two on site
         qo.push_back(Q(1,1));//a_up a_down

         for(int k = j - 1;k >= 0;--k){

            qo.push_back(Q(0,2));//a_down a_down
            qo.push_back(Q(1,1));//a_up a_down

         }

         for(int k = j - 1;k >= 0;--k){

            qo.push_back(Q(1,1));//a_up a_down
            qo.push_back(Q(2,0));//a_up a_up

         }

      }

      //a_0_up a_0_down
      qo.push_back(Q(1,1));

      for(int j = 2*i + 1;j < istates.size();++j)
         ostates.push_back(istates[j]);

      mpo[i].resize(Q::zero(),make_array(-mpo[i - 1].qshape(3),qp,-qp,qo));

      int column = 0;

      //identity
      mpo[i].insert(shape(0,0,0,0),Ip);
      mpo[i].insert(shape(0,1,1,0),Ip);
      mpo[i].insert(shape(0,2,2,0),Ip);
      mpo[i].insert(shape(0,3,3,0),Ip);

      //a_down
      mpo[i].insert(shape(0,0,2,1),Ip);
      mpo[i].insert(shape(0,1,3,1),Ip);

      //a_up (-1)^{n_down}
      mpo[i].insert(shape(0,0,1,2),Ip);
      mpo[i].insert(shape(0,2,3,2),Im);

      column = 3;

      //signs!
      for(int j = 1;j < 2*i + 1;++j){

         //fermion sign!
         mpo[i].insert(shape(j,0,0,column),Ip);
         mpo[i].insert(shape(j,1,1,column),Im);
         mpo[i].insert(shape(j,2,2,column),Im);
         mpo[i].insert(shape(j,3,3,column),Ip);

         column++;

      }

      //double remove
      mpo[i].insert(shape(0,0,3,column),Ip);
      column++;

      //single remove down (-1)^n_up
      for(int j = 1;j < 2*i + 1;++j){

         //fermion sign!
         mpo[i].insert(shape(j,0,2,column),Ip);
         mpo[i].insert(shape(j,1,3,column),Im);
         column++;

      }

      //single remove up
      for(int j = 1;j < 2*i + 1;++j){

         //fermion sign!
         mpo[i].insert(shape(j,0,1,column),Ip);
         mpo[i].insert(shape(j,2,3,column),Ip);
         ++column;

      }

      int row = 2*i + 1;

      for(int j = 0;j < i*(2*i - 1);++j){

         //finally identity
         mpo[i].insert(shape(row,0,0,column),Ip);
         mpo[i].insert(shape(row,1,1,column),Ip);
         mpo[i].insert(shape(row,2,2,column),Ip);
         mpo[i].insert(shape(row,3,3,column),Ip);
         ++row;
         ++column;

      }

      istates = ostates;

   }

   ostates.clear();

   //last occupied: i = no - 1
   qo.clear();

   qo.push_back(Q(1,1));//a_up a_down

   v.resize(4);
   v[0] = no - 1;
   v[1] = 1;
   v[2] = no - 1;
   v[3] = 0;

   ostates.push_back(v);

   for(int j = 1;j < no;++j){

      qo.push_back(Q(0,2));//a_down a_down

      v[0] = no - 1;
      v[1] = 1;
      v[2] = istates[2*j - 1][0];
      v[3] = istates[2*j - 1][1];

      ostates.push_back(v);

      qo.push_back(Q(1,1));//a_up a_down

      v[0] = no - 1;
      v[1] = 1;
      v[2] = istates[2*j][0];
      v[3] = istates[2*j][1];

      ostates.push_back(v);

   }

   for(int j = 1;j < no;++j){

      qo.push_back(Q(1,1));//a_down a_up

      v[0] = no - 1;
      v[1] = 0;
      v[2] = istates[2*j - 1][0];
      v[3] = istates[2*j - 1][1];

      ostates.push_back(v);

      qo.push_back(Q(2,0));//a_up a_up

      v[0] = no - 1;
      v[1] = 0;
      v[2] = istates[2*j][0];
      v[3] = istates[2*j][1];

      ostates.push_back(v);

   }

   //identity for the pairs from previous iteration
   for(int j = no - 2;j > 0;--j){

      //two on site
      qo.push_back(Q(1,1));//a_up a_down

      for(int k = j - 1;k >= 0;--k){

         qo.push_back(Q(0,2));//a_down a_down
         qo.push_back(Q(1,1));//a_down a_up

      }

      for(int k = j - 1;k >= 0;--k){

         qo.push_back(Q(1,1));//a_up a_down
         qo.push_back(Q(2,0));//a_up a_up

      } 

   }

   //a_0_up a_0_down
   qo.push_back(Q(1,1));

   for(int j = 2*no - 1;j < istates.size();++j)
      ostates.push_back(istates[j]);

   //start filling
   mpo[no - 1].resize(Q::zero(),make_array(-mpo[no - 2].qshape(3),qp,-qp,qo));

   int column = 0;

   //double remove
   mpo[no - 1].insert(shape(0,0,3,column),Ip);
   column++;

   //single remove down
   for(int j = 1;j < 2*no - 1;++j){

      //fermion sign!
      mpo[no - 1].insert(shape(j,0,2,column),Ip);
      mpo[no - 1].insert(shape(j,1,3,column),Im);
      column++;

   }

   //single remove up
   for(int j = 1;j < 2*no - 1;++j){

      mpo[no - 1].insert(shape(j,0,1,column),Ip);
      mpo[no - 1].insert(shape(j,2,3,column),Ip);
      ++column;

   }

   int row = 2*no - 1;

   for(int j = 0;j < (no - 1)*(2*(no - 1) - 1);++j){

      //finally identity
      mpo[no - 1].insert(shape(row,0,0,column),Ip);
      mpo[no - 1].insert(shape(row,1,1,column),Ip);
      mpo[no - 1].insert(shape(row,2,2,column),Ip);
      mpo[no - 1].insert(shape(row,3,3,column),Ip);

      ++row;
      ++column;

   }

   istates = ostates;

   ostates.clear();

   //first virtual
   qo.clear();

   qo.push_back(Q::zero());//first column completely closed 

   v.resize(1);
   v[0] = 1;
   ostates.push_back(v);

   //first all the virtuals
   for(int j = no + 1;j < L;++j){

      qo.push_back(Q(0,1));//go to 1 down removed

      v.resize(2);
      v[0] = j;
      v[1] = 1;
      ostates.push_back(v);

      qo.push_back(Q(1,0));//go to 1 up removed

      v[0] = j;
      v[1] = 0;
      ostates.push_back(v);

   }

   //the rest is identity
   for(int i = 0;i < istates.size();++i){

      ostates.push_back(istates[i]);
      qo.push_back(mpo[no - 1].qshape(3)[i]);

   }

   //start the filling
   mpo[no].resize(Q::zero(),make_array(-mpo[no - 1].qshape(3),qp,-qp,qo));

   //first the local term
   for(int j = 0;j < istates.size();++j){

      if(mpo[no].qshape(0)[j].gn_up() == -1 && mpo[no].qshape(0)[j].gn_down() == -1){//insert double create

         DArray<4> Tp(1,1,1,1);

         if(istates[j][1] == 1)//down spin is first
            Tp = t(istates[j][2],istates[j][0],0,0);
         else//up spin is first
            Tp = -t(istates[j][2],istates[j][0],0,0);

         //double create: first column
         mpo[no].insert(shape(j,3,0,0),Tp);

      }

   }

   //then the rest of the virtuals
   for(int j = no + 1;j < L;++j){

      int vjnd = j - no;//virtual index

      //first to one down removed
      for(int k = 0;k < istates.size();++k){//these are the rows, sum over pairs

         if(mpo[no].qshape(0)[k].gn_up() == -1 && mpo[no].qshape(0)[k].gn_down() == -1){//insert create up, so from (-1,-1) -> (0,1)

            DArray<4> Tp(1,1,1,1);
            DArray<4> Tm(1,1,1,1);

            if(istates[k][2] == istates[k][0]){//is it a double occupied site coming in?

               Tp = t(istates[k][2],istates[k][0],0,vjnd);
               Tm = -t(istates[k][2],istates[k][0],0,vjnd);

            }
            else{//two different sites coming in: up down or down up

               if(istates[k][1] == 1){// down up

                  Tp = t(istates[k][2],istates[k][0],0,vjnd);
                  Tm = -t(istates[k][2],istates[k][0],0,vjnd);

               }
               else{//up down

                  Tp = -t(istates[k][2],istates[k][0],vjnd,0);
                  Tm = t(istates[k][2],istates[k][0],vjnd,0);

               }

            }

            //create up: 
            mpo[no].insert(shape(k,1,0,2*vjnd-1),Tp);
            mpo[no].insert(shape(k,3,2,2*vjnd-1),Tm);

         }
         else if(mpo[no].qshape(0)[k].gn_up() == 0 && mpo[no].qshape(0)[k].gn_down() == -2){//insert create down: from (0,2) -> (0,-1)

            DArray<4> Tp(1,1,1,1);

            Tp = t(istates[k][2],istates[k][0],0,vjnd) - t(istates[k][2],istates[k][0],vjnd,0);

            //create down: 
            mpo[no].insert(shape(k,2,0,2*vjnd-1),Tp);
            mpo[no].insert(shape(k,3,1,2*vjnd-1),Tp);

         }

      }

      //then to one up removed
      for(int k = 0;k < istates.size();++k){

         if(mpo[no].qshape(0)[k].gn_up() == -1 && mpo[no].qshape(0)[k].gn_down() == -1){//insert create down, from (1,1) -> (-1,0)

            DArray<4> Tp(1,1,1,1);

            if(istates[k][2] == istates[k][0])//is it a double occupied site coming in?
               Tp = -t(istates[k][2],istates[k][0],0,vjnd);
            else{//two different sites coming in: up down or down up

               if(istates[k][1] == 1)// down up
                  Tp = -t(istates[k][2],istates[k][0],vjnd,0);
               else//up down
                  Tp = t(istates[k][2],istates[k][0],0,vjnd);

            }

            //create down: 
            mpo[no].insert(shape(k,2,0,2*vjnd),Tp);
            mpo[no].insert(shape(k,3,1,2*vjnd),Tp);

         }
         else if(mpo[no].qshape(0)[k].gn_up() == -2 && mpo[no].qshape(0)[k].gn_down() == 0){//insert create up : from (2,0) -> (-1,0)

            DArray<4> Tp(1,1,1,1);
            DArray<4> Tm(1,1,1,1);

            Tp = t(istates[k][2],istates[k][0],0,vjnd) - t(istates[k][2],istates[k][0],vjnd,0);
            Tm = t(istates[k][2],istates[k][0],vjnd,0) + t(istates[k][2],istates[k][0],0,vjnd);

            //create up:
            mpo[no].insert(shape(k,1,0,2*vjnd),Tp);
            mpo[no].insert(shape(k,3,2,2*vjnd),Tm);

         }

      }

   }

   //fill the remaining part unity
   for(int j = 0;j < istates.size();++j){

      mpo[no].insert(shape(j,0,0,j + 2*nv - 1),Ip);
      mpo[no].insert(shape(j,1,1,j + 2*nv - 1),Ip);
      mpo[no].insert(shape(j,2,2,j + 2*nv - 1),Ip);
      mpo[no].insert(shape(j,3,3,j + 2*nv - 1),Ip);

   }

   //this keeps track of all the pairs
   std::vector< std::vector<int> > pairs = istates;
   istates = ostates;

   //next all but the last virtuals:
   for(int i = no + 1;i < L - 1;++i){

      //current virtual index
      int vind = i - no;

      ostates.clear();

      qo.clear();

      qo.push_back(Q::zero());//first column completely closed 

      v.resize(1);
      v[0] = 1;
      ostates.push_back(v);

      //first all the remaining the virtuals
      for(int j = i + 1;j < L;++j){

         qo.push_back(Q(0,1));//go to 1 down removed

         v.resize(2);
         v[0] = j;
         v[1] = 1;
         ostates.push_back(v);

         qo.push_back(Q(1,0));//go to 1 down removed

         v[0] = j;
         v[1] = 0;
         ostates.push_back(v);

      }

      //then just the pairs
      for(int j = 0;j < pairs.size();++j){

         ostates.push_back(pairs[j]);
         qo.push_back(mpo[i - 1].qshape(3)[2*(nv - vind) + 1 + j]);

      }

      //start filling
      mpo[i].resize(Q::zero(),make_array(-mpo[i - 1].qshape(3),qp,-qp,qo));

      //first column is closed
      mpo[i].insert(shape(0,0,0,0),Ip);
      mpo[i].insert(shape(0,1,1,0),Ip);
      mpo[i].insert(shape(0,2,2,0),Ip);
      mpo[i].insert(shape(0,3,3,0),Ip);

      //create a down particle
      mpo[i].insert(shape(1,2,0,0),Ip);
      mpo[i].insert(shape(1,3,1,0),Im);

      //create an up particle
      mpo[i].insert(shape(2,1,0,0),Ip);
      mpo[i].insert(shape(2,3,2,0),Ip);

      row = 2*(nv - vind) + 1;

      //insert the closing pairs
      for(int j = 0;j < pairs.size();++j){

         if(mpo[i].qshape(0)[row].gn_up() == -1 && mpo[i].qshape(0)[row].gn_down() == -1){//insert double create

            DArray<4> Tp(1,1,1,1);

            if(pairs[j][1] == 1)//down spin is first
               Tp = t(pairs[j][2],pairs[j][0],vind,vind);
            else//up spin is first
               Tp = -t(pairs[j][2],pairs[j][0],vind,vind);

            //double create: first column
            mpo[i].insert(shape(row,3,0,0),Tp);

         }

         ++row;

      }

      //insert fermion signs
      for(int j = 3;j < 2*(nv - vind) + 1;++j){

         mpo[i].insert(shape(j,0,0,j-2),Ip);
         mpo[i].insert(shape(j,1,1,j-2),Im);
         mpo[i].insert(shape(j,2,2,j-2),Im);
         mpo[i].insert(shape(j,3,3,j-2),Ip);

      }

      //insert the t coefficients
      for(int j = i + 1;j < L;++j){

         int vjnd = j - no;
         int col = j - i;

         row = 2*(nv - vind) + 1;

         //first to one down removed
         for(int k = 0;k < pairs.size();++k){//these are the rows, sum over pairs

            if(mpo[i].qshape(0)[row].gn_up() == -1 && mpo[i].qshape(0)[row].gn_down() == -1){//insert create up, so from (-1,-1) -> (0,1)

               DArray<4> Tp(1,1,1,1);
               DArray<4> Tm(1,1,1,1);

               if(pairs[k][2] == pairs[k][0]){//is it a double occupied site coming in?

                  Tp = t(pairs[k][2],pairs[k][0],vind,vjnd);
                  Tm = -t(pairs[k][2],pairs[k][0],vind,vjnd);

               }
               else{//two different sites coming in: up down or down up

                  if(pairs[k][1] == 1){// down up

                     Tp = t(pairs[k][2],pairs[k][0],vind,vjnd);
                     Tm = -t(pairs[k][2],pairs[k][0],vind,vjnd);

                  }
                  else{//up down

                     Tp = -t(pairs[k][2],pairs[k][0],vjnd,vind);
                     Tm = t(pairs[k][2],pairs[k][0],vjnd,vind);

                  }

               }

               //create up: 
               mpo[i].insert(shape(row,1,0,2*col-1),Tp);
               mpo[i].insert(shape(row,3,2,2*col-1),Tm);

            }
            else if(mpo[i].qshape(0)[row].gn_up() == 0 && mpo[i].qshape(0)[row].gn_down() == -2){//insert create down: from (0,2) -> (0,-1)

               DArray<4> Tp(1,1,1,1);

               Tp = t(pairs[k][2],pairs[k][0],vind,vjnd) - t(pairs[k][2],pairs[k][0],vjnd,vind);

               //create down: 
               mpo[i].insert(shape(row,2,0,2*col-1),Tp);
               mpo[i].insert(shape(row,3,1,2*col-1),Tp);

            }

            ++row;

         }

         row = 2*(nv - vind) + 1;

         //then to one up removed
         for(int k = 0;k < pairs.size();++k){

            if(mpo[i].qshape(0)[row].gn_up() == -1 && mpo[i].qshape(0)[row].gn_down() == -1){//insert create down, from (1,1) -> (-1,0)

               DArray<4> Tp(1,1,1,1);

               if(pairs[k][2] == pairs[k][0])//is it a double occupied site coming in?
                  Tp = -t(pairs[k][2],pairs[k][0],vind,vjnd);
               else{//two different sites coming in: up down or down up

                  if(pairs[k][1] == 1)// down up
                     Tp = -t(pairs[k][2],pairs[k][0],vjnd,vind);
                  else//up down
                     Tp = t(pairs[k][2],pairs[k][0],vind,vjnd);

               }

               //create down: 
               mpo[i].insert(shape(row,2,0,2*col),Tp);
               mpo[i].insert(shape(row,3,1,2*col),Tp);

            }
            else if(mpo[i].qshape(0)[row].gn_up() == -2 && mpo[i].qshape(0)[row].gn_down() == 0){//insert create up : from (2,0) -> (-1,0)

               DArray<4> Tp(1,1,1,1);
               DArray<4> Tm(1,1,1,1);

               Tp = t(pairs[k][2],pairs[k][0],vind,vjnd) - t(pairs[k][2],pairs[k][0],vjnd,vind);
               Tm = t(pairs[k][2],pairs[k][0],vjnd,vind) + t(pairs[k][2],pairs[k][0],vind,vjnd);

               //create up:
               mpo[i].insert(shape(row,1,0,2*col),Tp);
               mpo[i].insert(shape(row,3,2,2*col),Tm);

            }

            ++row;

         }

      }

      //fill the lower right block with unity
      for(int j = 2*(nv - vind) + 1;j < istates.size();++j){

         mpo[i].insert(shape(j,0,0,j - 2),Ip);
         mpo[i].insert(shape(j,1,1,j - 2),Ip);
         mpo[i].insert(shape(j,2,2,j - 2),Ip);
         mpo[i].insert(shape(j,3,3,j - 2),Ip);

      }

      istates = ostates;

   }

   //last virtual L - 1, closing everything down
   ostates.clear();
   qo.clear();

   qo.push_back(Q::zero());//first column completely closed 

   v.resize(1);
   v[0] = 1;
   ostates.push_back(v);

   mpo[L - 1].resize(Q::zero(),make_array(-mpo[L - 2].qshape(3),qp,-qp,qo));

   //first column is closed
   mpo[L - 1].insert(shape(0,0,0,0),Ip);
   mpo[L - 1].insert(shape(0,1,1,0),Ip);
   mpo[L - 1].insert(shape(0,2,2,0),Ip);
   mpo[L - 1].insert(shape(0,3,3,0),Ip);

   //create a down particle
   mpo[L - 1].insert(shape(1,2,0,0),Ip);
   mpo[L - 1].insert(shape(1,3,1,0),Im);

   //create an up particle
   mpo[L - 1].insert(shape(2,1,0,0),Ip);
   mpo[L - 1].insert(shape(2,3,2,0),Ip);

   row = 3;

   //insert the closing pairs
   for(int j = 0;j < pairs.size();++j){

      if(mpo[L - 1].qshape(0)[row].gn_up() == -1 && mpo[L - 1].qshape(0)[row].gn_down() == -1){//insert double create

         DArray<4> Tp(1,1,1,1);

         if(pairs[j][1] == 1)//down spin is first
            Tp = t(pairs[j][2],pairs[j][0],nv - 1,nv - 1);
         else//up spin is first
            Tp = -t(pairs[j][2],pairs[j][0],nv - 1,nv - 1);

         //double create: first column
         mpo[L - 1].insert(shape(row,3,0,0),Tp);

      }

      ++row;

   }

   //merge everything together
   TVector<Qshapes<Q>,1> qmerge;
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
*/
   return mpo;

}

template void physical<Quantum>(Qshapes<Quantum> &);
template MPO<Quantum> creator<Quantum>(int L,int site,int spin);
template MPO<Quantum> annihilator(int L,int site,int spin);
template MPO<Quantum> n_loc(int L,int site);
template MPO<Quantum> N_tot(int L);
template MPO<Quantum> n_up_tot(int L);
template MPO<Quantum> n_down_tot(int L);
template MPO<Quantum> hubbard(int L,double U);
template MPO<Quantum> T1(const DArray<2> &);
template MPO<Quantum> T2(const DArray<4> &);
template MPO<Quantum> one_body(const DArray<2> &);
