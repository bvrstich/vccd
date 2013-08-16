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

   qp.push_back(Q(0,0));
   qp.push_back(Q(0,1));
   qp.push_back(Q(1,0));
   qp.push_back(Q(1,1));

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

template void physical<Quantum>(Qshapes<Quantum> &);
template MPO<Quantum> creator<Quantum>(int L,int site,int spin);
template MPO<Quantum> annihilator(int L,int site,int spin);
template MPO<Quantum> n_loc(int L,int site);
template MPO<Quantum> N_tot(int L);
template MPO<Quantum> n_up_tot(int L);
template MPO<Quantum> n_down_tot(int L);
template MPO<Quantum> hubbard(int L,double U);
template MPO<Quantum> T1(const DArray<2> &);
