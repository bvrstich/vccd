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
 * MPO representing creation operator of up spin paritcle on site i
 */
template<class Q>
MPO<Q> crea_up(int L,int i){

   MPO<Q> mpo(L);

   Qshapes<Q> qp;
   physical(qp);

   Qshapes<Q> qz; // 0 quantum number
   qz.push_back(Q::zero());

   //signs before the up
   for(int site = 0;site < i;++site){

      mpo[site].resize(Q::zero(),make_array(qz,qp,-qp,qz));
      insert_sign(mpo[site],0,0);

   }

   //create up spin on site i
   Qshapes<Q> qo; 
   qo.push_back(Q(-1,0));

   mpo[i].resize(Q::zero(),make_array(qz,qp,-qp,qo));
   insert_crea_up(mpo[i],0,0,1.0);

   //identitity after the operator
   for(int site = i + 1;site < L;++site){

      mpo[site].resize(Q::zero(),make_array(-qo,qp,-qp,qo));
      insert_id(mpo[site],0,0);

   }

   return mpo;

}

/**
 * MPO representing creation operator of a down spin paritcle on site i
 */
template<class Q>
MPO<Q> crea_down(int L,int i){

   MPO<Q> mpo(L);

   Qshapes<Q> qp;
   physical(qp);

   Qshapes<Q> qz; // 0 quantum number
   qz.push_back(Q::zero());

   //signs before the up
   for(int site = 0;site < i;++site){

      mpo[site].resize(Q::zero(),make_array(qz,qp,-qp,qz));
      insert_sign(mpo[site],0,0);

   }

   //create down spin on site i
   Qshapes<Q> qo; 
   qo.push_back(Q(0,-1));

   mpo[i].resize(Q::zero(),make_array(qz,qp,-qp,qo));
   insert_crea_down_s(mpo[i],0,0,1.0);

   //identitity after the operator
   for(int site = i + 1;site < L;++site){

      mpo[site].resize(Q::zero(),make_array(-qo,qp,-qp,qo));
      insert_id(mpo[site],0,0);

   }

   return mpo;

}

/**
 * MPO representing annihilation operator of up spin paritcle on site i
 */
template<class Q>
MPO<Q> anni_up(int L,int i){

   MPO<Q> mpo(L);

   Qshapes<Q> qp;
   physical(qp);

   Qshapes<Q> qz; // 0 quantum number
   qz.push_back(Q::zero());

   //signs before the up
   for(int site = 0;site < i;++site){

      mpo[site].resize(Q::zero(),make_array(qz,qp,-qp,qz));
      insert_sign(mpo[site],0,0);

   }

   //anni up spin on site i
   Qshapes<Q> qo; 
   qo.push_back(Q(1,0));

   mpo[i].resize(Q::zero(),make_array(qz,qp,-qp,qo));
   insert_anni_up(mpo[i],0,0,1.0);

   //identitity after the operator
   for(int site = i + 1;site < L;++site){

      mpo[site].resize(Q::zero(),make_array(-qo,qp,-qp,qo));
      insert_id(mpo[site],0,0);

   }

   return mpo;

}

/**
 * MPO representing annihilation operator of a down spin paritcle on site i
 */
template<class Q>
MPO<Q> anni_down(int L,int i){

   MPO<Q> mpo(L);

   Qshapes<Q> qp;
   physical(qp);

   Qshapes<Q> qz; // 0 quantum number
   qz.push_back(Q::zero());

   //signs before the down
   for(int site = 0;site < i;++site){

      mpo[site].resize(Q::zero(),make_array(qz,qp,-qp,qz));
      insert_sign(mpo[site],0,0);

   }

   //annihilate down spin on site i
   Qshapes<Q> qo; 
   qo.push_back(Q(0,1));

   mpo[i].resize(Q::zero(),make_array(qz,qp,-qp,qo));
   insert_anni_down_s(mpo[i],0,0,1.0);

   //identitity after the operator
   for(int site = i + 1;site < L;++site){

      mpo[site].resize(Q::zero(),make_array(-qo,qp,-qp,qo));
      insert_id(mpo[site],0,0);

   }

   return mpo;

}

/**
 * MPO representing the creation operator of a pair of up spin paritcle on site i and j
 */
template<class Q>
MPO<Q> crea_up_crea_up(int L,int i,int j,double val){

   MPO<Q> mpo(L);

   Qshapes<Q> qp;
   physical(qp);

   Qshapes<Q> qz; // 0 quantum number
   qz.push_back(Q::zero());

   if(i == j){//just zero everywhere

      for(int site = 0;site < L;++site){

         mpo[site].resize(Q::zero(),make_array(qz,qp,-qp,qz));
         insert_zero(mpo[site],0,0);

      }

   }
   else{

      int sign = 1;

      if(i > j){

         sign *= -1;
         int tmp = i;
         i = j;
         j = tmp;

      }

      //id before the first up
      for(int site = 0;site < i;++site){

         mpo[site].resize(Q::zero(),make_array(qz,qp,-qp,qz));
         insert_id(mpo[site],0,0);

      }

      //create up spin on site i
      Qshapes<Q> qo; 
      qo.push_back(Q(-1,0));

      mpo[i].resize(Q::zero(),make_array(qz,qp,-qp,qo));
      insert_crea_up_s(mpo[i],0,0,val);

      //signs between
      for(int site = i + 1;site < j;++site){

         mpo[site].resize(Q::zero(),make_array(-qo,qp,-qp,qo));
         insert_sign(mpo[site],0,0);

      }

      //create up spin on site j
      Qshapes<Q> qi = qo; 

      qo.clear();
      qo.push_back(Q(-2,0));

      mpo[j].resize(Q::zero(),make_array(-qi,qp,-qp,qo));
      insert_crea_up(mpo[j],0,0,sign);

      //id after j
      for(int site = j + 1;site < L;++site){

         mpo[site].resize(Q::zero(),make_array(-qo,qp,-qp,qo));
         insert_id(mpo[site],0,0);

      }

   }

   return mpo;

}

/**
 * MPO representing the creation operator of a pair of up spin paritcle on site i and j
 */
template<class Q>
MPO<Q> crea_up_crea_down(int L,int i,int j,double val){

   MPO<Q> mpo(L);

   Qshapes<Q> qp;
   physical(qp);

   Qshapes<Q> qz; // 0 quantum number
   qz.push_back(Q::zero());

   if(i < j){

      //id before the first up
      for(int site = 0;site < i;++site){

         mpo[site].resize(Q::zero(),make_array(qz,qp,-qp,qz));
         insert_id(mpo[site],0,0);

      }

      //create up spin on site i
      Qshapes<Q> qo; 
      qo.push_back(Q(-1,0));

      mpo[i].resize(Q::zero(),make_array(qz,qp,-qp,qo));
      insert_crea_up_s(mpo[i],0,0,val);

      //signs between i and j
      for(int site = i + 1;site < j;++site){

         mpo[site].resize(Q::zero(),make_array(-qo,qp,-qp,qo));
         insert_sign(mpo[site],0,0);

      }

      //create down spin on site j
      Qshapes<Q> qi = qo; 

      qo.clear();
      qo.push_back(Q(-1,-1));

      mpo[j].resize(Q::zero(),make_array(-qi,qp,-qp,qo));
      insert_crea_down_s(mpo[j],0,0,1.0);

      //id after j
      for(int site = j + 1;site < L;++site){

         mpo[site].resize(Q::zero(),make_array(-qo,qp,-qp,qo));
         insert_id(mpo[site],0,0);

      }

   }
   else if(i == j){

      //id before the operator
      for(int site = 0;site < i;++site){

         mpo[site].resize(Q::zero(),make_array(qz,qp,-qp,qz));
         insert_id(mpo[site],0,0);

      }

      //create up and down on site i
      Qshapes<Q> qo; 
      qo.push_back(Q(-1,-1));

      mpo[i].resize(Q::zero(),make_array(qz,qp,-qp,qo));
      insert_crea_up_crea_down(mpo[i],0,0,val);

      //id after
      for(int site = i + 1;site < L;++site){

         mpo[site].resize(Q::zero(),make_array(-qo,qp,-qp,qo));
         insert_id(mpo[site],0,0);

      }

   }
   else{//i > j

      //id before the first down
      for(int site = 0;site < j;++site){

         mpo[site].resize(Q::zero(),make_array(qz,qp,-qp,qz));
         insert_id(mpo[site],0,0);

      }

      //create down spin on site j
      Qshapes<Q> qo; 
      qo.push_back(Q(0,-1));

      mpo[j].resize(Q::zero(),make_array(qz,qp,-qp,qo));
      insert_crea_down(mpo[j],0,0,val);

      //signs between j and i
      for(int site = j + 1;site < i;++site){

         mpo[site].resize(Q::zero(),make_array(-qo,qp,-qp,qo));
         insert_sign(mpo[site],0,0);

      }

      //create up spin on site i
      Qshapes<Q> qi = qo; 

      qo.clear();
      qo.push_back(Q(-1,-1));

      mpo[i].resize(Q::zero(),make_array(-qi,qp,-qp,qo));
      insert_crea_up(mpo[i],0,0,-1.0);

      //id after i
      for(int site = i + 1;site < L;++site){

         mpo[site].resize(Q::zero(),make_array(-qo,qp,-qp,qo));
         insert_id(mpo[site],0,0);

      }

   }

   return mpo;

}

/**
 * MPO representing the creation operator of a pair of spin up and down paritcle on site i and j
 * i = down, j = up
 */
template<class Q>
MPO<Q> crea_down_crea_up(int L,int i,int j,double val){

   MPO<Q> mpo(L);

   Qshapes<Q> qp;
   physical(qp);

   Qshapes<Q> qz; // 0 quantum number
   qz.push_back(Q::zero());

   if(i < j){

      //id before the first up
      for(int site = 0;site < i;++site){

         mpo[site].resize(Q::zero(),make_array(qz,qp,-qp,qz));
         insert_id(mpo[site],0,0);

      }

      //create down spin on site i
      Qshapes<Q> qo; 
      qo.push_back(Q(0,-1));

      mpo[i].resize(Q::zero(),make_array(qz,qp,-qp,qo));
      insert_crea_down(mpo[i],0,0,val);

      //signs between i and j
      for(int site = i + 1;site < j;++site){

         mpo[site].resize(Q::zero(),make_array(-qo,qp,-qp,qo));
         insert_sign(mpo[site],0,0);

      }

      //create up spin on site j
      Qshapes<Q> qi = qo; 

      qo.clear();
      qo.push_back(Q(-1,-1));

      mpo[j].resize(Q::zero(),make_array(-qi,qp,-qp,qo));
      insert_crea_up(mpo[j],0,0,1.0);

      //id after j
      for(int site = j + 1;site < L;++site){

         mpo[site].resize(Q::zero(),make_array(-qo,qp,-qp,qo));
         insert_id(mpo[site],0,0);

      }

   }
   else if(i == j){

      //id before the operator
      for(int site = 0;site < i;++site){

         mpo[site].resize(Q::zero(),make_array(qz,qp,-qp,qz));
         insert_id(mpo[site],0,0);

      }

      //create up and down on site i
      Qshapes<Q> qo; 
      qo.push_back(Q(-1,-1));

      mpo[i].resize(Q::zero(),make_array(qz,qp,-qp,qo));
      insert_crea_up_crea_down(mpo[i],0,0,-val);

      //id after
      for(int site = i + 1;site < L;++site){

         mpo[site].resize(Q::zero(),make_array(-qo,qp,-qp,qo));
         insert_id(mpo[site],0,0);

      }

   }
   else{//i > j

      //id before the first down
      for(int site = 0;site < j;++site){

         mpo[site].resize(Q::zero(),make_array(qz,qp,-qp,qz));
         insert_id(mpo[site],0,0);

      }

      //create up spin on site j
      Qshapes<Q> qo; 
      qo.push_back(Q(-1,0));

      mpo[j].resize(Q::zero(),make_array(qz,qp,-qp,qo));
      insert_crea_up_s(mpo[j],0,0,val);

      //signs between j and i
      for(int site = j + 1;site < i;++site){

         mpo[site].resize(Q::zero(),make_array(-qo,qp,-qp,qo));
         insert_sign(mpo[site],0,0);

      }

      //create down spin on site i
      Qshapes<Q> qi = qo; 

      qo.clear();
      qo.push_back(Q(-1,-1));

      mpo[i].resize(Q::zero(),make_array(-qi,qp,-qp,qo));
      insert_crea_down_s(mpo[i],0,0,-1.0);

      //id after i
      for(int site = i + 1;site < L;++site){

         mpo[site].resize(Q::zero(),make_array(-qo,qp,-qp,qo));
         insert_id(mpo[site],0,0);

      }

   }

   return mpo;

}

/**
 * MPO representing the creation operator of a pair of up spin paritcle on site i and j
 */
template<class Q>
MPO<Q> crea_down_crea_down(int L,int i,int j,double val){

   MPO<Q> mpo(L);

   Qshapes<Q> qp;
   physical(qp);

   Qshapes<Q> qz; // 0 quantum number
   qz.push_back(Q::zero());

   if(i == j){//just zero everywhere

      for(int site = 0;site < L;++site){

         mpo[site].resize(Q::zero(),make_array(qz,qp,-qp,qz));
         insert_zero(mpo[site],0,0);

      }

   }
   else{

      int sign = 1;

      if(i > j){

         sign *= -1;
         int tmp = i;
         i = j;
         j = tmp;

      }

      //id before the first down
      for(int site = 0;site < i;++site){

         mpo[site].resize(Q::zero(),make_array(qz,qp,-qp,qz));
         insert_id(mpo[site],0,0);

      }

      //create up spin on site i
      Qshapes<Q> qo; 
      qo.push_back(Q(0,-1));

      mpo[i].resize(Q::zero(),make_array(qz,qp,-qp,qo));
      insert_crea_down(mpo[i],0,0,val);

      //signs between
      for(int site = i + 1;site < j;++site){

         mpo[site].resize(Q::zero(),make_array(-qo,qp,-qp,qo));
         insert_sign(mpo[site],0,0);

      }

      //create up spin on site j
      Qshapes<Q> qi = qo; 

      qo.clear();
      qo.push_back(Q(0,-2));

      mpo[j].resize(Q::zero(),make_array(-qi,qp,-qp,qo));
      insert_crea_down_s(mpo[j],0,0,sign);

      //id after j
      for(int site = j + 1;site < L;++site){

         mpo[site].resize(Q::zero(),make_array(-qo,qp,-qp,qo));
         insert_id(mpo[site],0,0);

      }

   }

   return mpo;

}

/**
 * MPO representing the creation operator of a pair of up spin paritcle on site i and j
 */
template<class Q>
MPO<Q> anni_up_anni_up(int L,int i,int j,double val){

   MPO<Q> mpo(L);

   Qshapes<Q> qp;
   physical(qp);

   Qshapes<Q> qz; // 0 quantum number
   qz.push_back(Q::zero());

   if(i == j){//just zero everywhere

      for(int site = 0;site < L;++site){

         mpo[site].resize(Q::zero(),make_array(qz,qp,-qp,qz));
         insert_zero(mpo[site],0,0);

      }

   }
   else{

      int sign = -1;

      if(i > j){

         sign *= -1;
         int tmp = i;
         i = j;
         j = tmp;

      }

      //id before the first up
      for(int site = 0;site < i;++site){

         mpo[site].resize(Q::zero(),make_array(qz,qp,-qp,qz));
         insert_id(mpo[site],0,0);

      }

      //create up spin on site i
      Qshapes<Q> qo; 
      qo.push_back(Q(1,0));

      mpo[i].resize(Q::zero(),make_array(qz,qp,-qp,qo));
      insert_anni_up_s(mpo[i],0,0,val);

      //signs between
      for(int site = i + 1;site < j;++site){

         mpo[site].resize(Q::zero(),make_array(-qo,qp,-qp,qo));
         insert_sign(mpo[site],0,0);

      }

      //anni up spin on site j
      Qshapes<Q> qi = qo; 

      qo.clear();
      qo.push_back(Q(2,0));

      mpo[j].resize(Q::zero(),make_array(-qi,qp,-qp,qo));
      insert_anni_up(mpo[j],0,0,sign);

      //id after j
      for(int site = j + 1;site < L;++site){

         mpo[site].resize(Q::zero(),make_array(-qo,qp,-qp,qo));
         insert_id(mpo[site],0,0);

      }

   }

   return mpo;

}

/**
 * MPO representing the annihilation operator of a down and a up spin
 */
template<class Q>
MPO<Q> anni_down_anni_up(int L,int i,int j,double val){

   MPO<Q> mpo(L);

   Qshapes<Q> qp;
   physical(qp);

   Qshapes<Q> qz; // 0 quantum number
   qz.push_back(Q::zero());

   if(i < j){

      //id before the first up
      for(int site = 0;site < i;++site){

         mpo[site].resize(Q::zero(),make_array(qz,qp,-qp,qz));
         insert_id(mpo[site],0,0);

      }

      //anni down spin on site i
      Qshapes<Q> qo; 
      qo.push_back(Q(0,1));

      mpo[i].resize(Q::zero(),make_array(qz,qp,-qp,qo));
      insert_anni_down(mpo[i],0,0,val);

      //signs between i and j
      for(int site = i + 1;site < j;++site){

         mpo[site].resize(Q::zero(),make_array(-qo,qp,-qp,qo));
         insert_sign(mpo[site],0,0);

      }

      //anni up spin on site j
      Qshapes<Q> qi = qo; 

      qo.clear();
      qo.push_back(Q(1,1));

      mpo[j].resize(Q::zero(),make_array(-qi,qp,-qp,qo));
      insert_anni_up(mpo[j],0,0,-1.0);

      //id after j
      for(int site = j + 1;site < L;++site){

         mpo[site].resize(Q::zero(),make_array(-qo,qp,-qp,qo));
         insert_id(mpo[site],0,0);

      }

   }
   else if(i == j){

      //id before the operator
      for(int site = 0;site < i;++site){

         mpo[site].resize(Q::zero(),make_array(qz,qp,-qp,qz));
         insert_id(mpo[site],0,0);

      }

      //anni down and up on site i
      Qshapes<Q> qo; 
      qo.push_back(Q(1,1));

      //extra minus sign!
      mpo[i].resize(Q::zero(),make_array(qz,qp,-qp,qo));
      insert_anni_down_anni_up(mpo[i],0,0,-val);

      //id after
      for(int site = i + 1;site < L;++site){

         mpo[site].resize(Q::zero(),make_array(-qo,qp,-qp,qo));
         insert_id(mpo[site],0,0);

      }

   }
   else{//i > j

      //id before the first up
      for(int site = 0;site < j;++site){

         mpo[site].resize(Q::zero(),make_array(qz,qp,-qp,qz));
         insert_id(mpo[site],0,0);

      }

      //anni up spin on site j
      Qshapes<Q> qo; 
      qo.push_back(Q(1,0));

      mpo[j].resize(Q::zero(),make_array(qz,qp,-qp,qo));
      insert_anni_up_s(mpo[j],0,0,val);

      //signs between j and i
      for(int site = j + 1;site < i;++site){

         mpo[site].resize(Q::zero(),make_array(-qo,qp,-qp,qo));
         insert_sign(mpo[site],0,0);

      }

      //anni down spin on site i
      Qshapes<Q> qi = qo; 

      qo.clear();
      qo.push_back(Q(1,1));

      mpo[i].resize(Q::zero(),make_array(-qi,qp,-qp,qo));
      insert_anni_down_s(mpo[i],0,0,1.0);

      //id after i
      for(int site = i + 1;site < L;++site){

         mpo[site].resize(Q::zero(),make_array(-qo,qp,-qp,qo));
         insert_id(mpo[site],0,0);

      }

   }

   return mpo;


}

/**
 * MPO representing the annihilation operator of an up and a down spin
 */
template<class Q>
MPO<Q> anni_up_anni_down(int L,int i,int j,double val){

   MPO<Q> mpo(L);

   Qshapes<Q> qp;
   physical(qp);

   Qshapes<Q> qz; // 0 quantum number
   qz.push_back(Q::zero());

   if(i < j){

      //id before the first up
      for(int site = 0;site < i;++site){

         mpo[site].resize(Q::zero(),make_array(qz,qp,-qp,qz));
         insert_id(mpo[site],0,0);

      }

      //anni up spin on site i
      Qshapes<Q> qo; 
      qo.push_back(Q(1,0));

      mpo[i].resize(Q::zero(),make_array(qz,qp,-qp,qo));
      insert_anni_up_s(mpo[i],0,0,val);

      //signs between i and j
      for(int site = i + 1;site < j;++site){

         mpo[site].resize(Q::zero(),make_array(-qo,qp,-qp,qo));
         insert_sign(mpo[site],0,0);

      }

      //anni up spin on site j
      Qshapes<Q> qi = qo; 

      qo.clear();
      qo.push_back(Q(1,1));

      mpo[j].resize(Q::zero(),make_array(-qi,qp,-qp,qo));
      insert_anni_down_s(mpo[j],0,0,-1.0);

      //id after j
      for(int site = j + 1;site < L;++site){

         mpo[site].resize(Q::zero(),make_array(-qo,qp,-qp,qo));
         insert_id(mpo[site],0,0);

      }

   }
   else if(i == j){

      //id before the operator
      for(int site = 0;site < i;++site){

         mpo[site].resize(Q::zero(),make_array(qz,qp,-qp,qz));
         insert_id(mpo[site],0,0);

      }

      //anni down and up on site i
      Qshapes<Q> qo; 
      qo.push_back(Q(1,1));

      mpo[i].resize(Q::zero(),make_array(qz,qp,-qp,qo));
      insert_anni_down_anni_up(mpo[i],0,0,val);

      //id after
      for(int site = i + 1;site < L;++site){

         mpo[site].resize(Q::zero(),make_array(-qo,qp,-qp,qo));
         insert_id(mpo[site],0,0);

      }

   }
   else{//i > j

      //id before the first up
      for(int site = 0;site < j;++site){

         mpo[site].resize(Q::zero(),make_array(qz,qp,-qp,qz));
         insert_id(mpo[site],0,0);

      }

      //anni down spin on site j
      Qshapes<Q> qo; 
      qo.push_back(Q(0,1));

      mpo[j].resize(Q::zero(),make_array(qz,qp,-qp,qo));
      insert_anni_down(mpo[j],0,0,val);

      //signs between j and i
      for(int site = j + 1;site < i;++site){

         mpo[site].resize(Q::zero(),make_array(-qo,qp,-qp,qo));
         insert_sign(mpo[site],0,0);

      }

      //anni up spin on site i
      Qshapes<Q> qi = qo; 

      qo.clear();
      qo.push_back(Q(1,1));

      mpo[i].resize(Q::zero(),make_array(-qi,qp,-qp,qo));
      insert_anni_up(mpo[i],0,0,1.0);

      //id after i
      for(int site = i + 1;site < L;++site){

         mpo[site].resize(Q::zero(),make_array(-qo,qp,-qp,qo));
         insert_id(mpo[site],0,0);

      }

   }

   return mpo;


}

/**
 * MPO representing the annihilation operator of a pair of down spin paritcle on site i and j
 */
template<class Q>
MPO<Q> anni_down_anni_down(int L,int i,int j,double val){

   MPO<Q> mpo(L);

   Qshapes<Q> qp;
   physical(qp);

   Qshapes<Q> qz; // 0 quantum number
   qz.push_back(Q::zero());

   if(i == j){//just zero everywhere

      for(int site = 0;site < L;++site){

         mpo[site].resize(Q::zero(),make_array(qz,qp,-qp,qz));
         insert_zero(mpo[site],0,0);

      }

   }
   else{

      int sign = -1;

      if(i > j){

         sign *= -1;
         int tmp = i;
         i = j;
         j = tmp;

      }

      //id before the first up
      for(int site = 0;site < i;++site){

         mpo[site].resize(Q::zero(),make_array(qz,qp,-qp,qz));
         insert_id(mpo[site],0,0);

      }

      //anni down spin on site i
      Qshapes<Q> qo; 
      qo.push_back(Q(0,1));

      mpo[i].resize(Q::zero(),make_array(qz,qp,-qp,qo));
      insert_anni_down(mpo[i],0,0,val);

      //signs between
      for(int site = i + 1;site < j;++site){

         mpo[site].resize(Q::zero(),make_array(-qo,qp,-qp,qo));
         insert_sign(mpo[site],0,0);

      }

      //anni up spin on site j
      Qshapes<Q> qi = qo; 

      qo.clear();
      qo.push_back(Q(0,2));

      mpo[j].resize(Q::zero(),make_array(-qi,qp,-qp,qo));
      insert_anni_down_s(mpo[j],0,0,sign);

      //id after j
      for(int site = j + 1;site < L;++site){

         mpo[site].resize(Q::zero(),make_array(-qo,qp,-qp,qo));
         insert_id(mpo[site],0,0);

      }

   }

   return mpo;

}

/**
 * general two-particle interaction term: represent as MPO
 */
template<class Q>
MPO<Q> tpint(int L,int i,int j,int k,int l,double V){

   MPO<Q> mpo;

   if(i != j && k != l){

      //up up up up
      MPO<Q> c = crea_up_crea_up<Q>(L,i,j,1.0);
      MPO<Q> a = anni_up_anni_up<Q>(L,l,k,V);

      mpo = c * a;

      //down down down down
      c = crea_down_crea_down<Q>(L,i,j,1.0);
      a = anni_down_anni_down<Q>(L,l,k,V);

      gemm(1.0,c,a,1.0,mpo);

      //down up up down
      c = crea_down_crea_up<Q>(L,i,j,1.0);
      a = anni_up_anni_down<Q>(L,l,k,V);

      gemm(1.0,c,a,1.0,mpo);

      //down up up down
      c = crea_up_crea_down<Q>(L,i,j,1.0);
      a = anni_down_anni_up<Q>(L,l,k,V);

      gemm(1.0,c,a,1.0,mpo);

      return mpo;

   }
   else{

      //down up up down
      MPO<Q> c = crea_down_crea_up<Q>(L,i,j,1.0);
      MPO<Q> a = anni_up_anni_down<Q>(L,l,k,V);

      mpo = c*a;

      //down up up down
      c = crea_up_crea_down<Q>(L,i,j,1.0);
      a = anni_down_anni_up<Q>(L,l,k,V);

      gemm(1.0,c,a,1.0,mpo);

      return mpo;

   }

}

/**
 * elementary double excitation operator:E^i_k E^j_l
 */
template<class Q>
MPO<Q> E(int L,int i,int j,int k,int l,double t){

   MPO<Q> Eik = E<Q>(L,i,k,1.0);
   MPO<Q> Ejl = E<Q>(L,j,l,t);

   return Eik*Ejl;

}

/**
 * elementary double excitation operator:E^i_k E^j_l
 */
template<class Q>
MPO<Q> T2_test(const DArray<4> &t){

   int no = t.shape(0);
   int nv = t.shape(2);

   int L = no + nv;

   MPO<Q> tmp = E<Q>(L,no,no,0,0,t(0,0,0,0));

   for(int i = 0;i < no;++i)
      for(int j = 0;j < no;++j)
         for(int a = 0;a < nv;++a)
            for(int b = 0;b < nv;++b){

               if( (i + j + a + b) != 0 )
                  mpsxx::axpy(1.0,E<Q>(L,a + no,b + no,i,j,t(i,j,a,b)),tmp);

            }

   return tmp;

}

/**
 * elementary excitation operator:E^i_j: \sum_s a^+_{i,s} a_{js}
 */
template<class Q>
MPO<Q> E(int L,int i,int j,double t){

   MPO<Q> mpo(L);

   Qshapes<Q> qp;
   physical(qp);

   Qshapes<Q> qz; // 0 quantum number
   qz.push_back(Q::zero());

   if(i < j){//first creator

      //left of i
      for(int site = 0;site < i;++site){

         mpo[site].resize(Q::zero(),make_array(qz,qp,-qp,qz));
         insert_id(mpo[site],0,0);

      }

      //site == i: insert creator

      Qshapes<Q> qo; // 0 quantum number
      qo.push_back(Q(-1,0));
      qo.push_back(Q(0,-1));

      mpo[i].resize(Q::zero(),make_array(qz,qp,-qp,qo));

      insert_crea_up_s(mpo[i],0,0,1.0);
      insert_crea_down(mpo[i],0,1,1.0);

      //signs in between
      for(int site = i + 1;site < j;++site){

         mpo[site].resize(Q::zero(),make_array(-qo,qp,-qp,qo));

         insert_sign(mpo[site],0,0);
         insert_sign(mpo[site],1,1);

      }

      //close down on site j
      mpo[j].resize(Q::zero(),make_array(-qo,qp,-qp,qz));

      insert_anni_up(mpo[j],0,0,t);
      insert_anni_down_s(mpo[j],1,0,t);

      //id on the rest of lattice
      for(int site = j + 1;site < L;++site){

         mpo[site].resize(Q::zero(),make_array(qz,qp,-qp,qz));
         insert_id(mpo[site],0,0);

      }

   }
   else if(i == j){//just local term

      //left of i
      for(int site = 0;site < i;++site){

         mpo[site].resize(Q::zero(),make_array(qz,qp,-qp,qz));
         insert_id(mpo[site],0,0);

      }

      //on i
      mpo[i].resize(Q::zero(),make_array(qz,qp,-qp,qz));
      insert_local_ob(mpo[i],0,0,t);

      //right of i
      for(int site = i + 1;site < L;++site){

         mpo[site].resize(Q::zero(),make_array(qz,qp,-qp,qz));
         insert_id(mpo[site],0,0);

      }

   }
   else{//i > j

      //left of j
      for(int site = 0;site < j;++site){

         mpo[site].resize(Q::zero(),make_array(qz,qp,-qp,qz));
         insert_id(mpo[site],0,0);

      }

      //site == j: insert annihilator

      Qshapes<Q> qo; // 0 quantum number
      qo.push_back(Q(1,0));
      qo.push_back(Q(0,1));

      mpo[j].resize(Q::zero(),make_array(qz,qp,-qp,qo));

      insert_anni_up_s(mpo[j],0,0,1.0);
      insert_anni_down(mpo[j],0,1,1.0);

      //signs in between
      for(int site = j + 1;site < i;++site){

         mpo[site].resize(Q::zero(),make_array(-qo,qp,-qp,qo));

         insert_sign(mpo[site],0,0);
         insert_sign(mpo[site],1,1);

      }

      //close down on site i
      mpo[i].resize(Q::zero(),make_array(-qo,qp,-qp,qz));

      insert_crea_up(mpo[i],0,0,t);
      insert_crea_down_s(mpo[i],1,0,t);

      //id on the rest of lattice
      for(int site = i + 1;site < L;++site){

         mpo[site].resize(Q::zero(),make_array(qz,qp,-qp,qz));
         insert_id(mpo[site],0,0);

      }

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
 * @return MPO object of length L containing the T1 operator with coefficients passed through the DArray<2> object t
 */
template<class Q>
MPO<Q> T1(const DArray<2> &t,bool merge){

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
   qo.push_back(Q(0,1));//a_down
   qo.push_back(Q(1,0));//a_up

   std::vector< Ostate > ostates;
   Ostate state;

   state.push_id();
   ostates.push_back(state);
   state.clear();

   state.push_anni_down(0);
   ostates.push_back(state);
   state.clear();

   state.push_anni_up(0);
   ostates.push_back(state);
   state.clear();

   mpo[0].resize(Q::zero(),make_array(qz,qp,-qp,qo));

   insert_id(mpo[0],0,0);
   insert_anni_down(mpo[0],0,1,1.0);
   insert_anni_up_s(mpo[0],0,2,1.0);

   std::vector< Ostate > istates;
   istates = ostates;

   Qshapes<Quantum> qi = qo;

   for(int i = 1;i < no - 1;++i){

      ostates.clear();
      qo.clear();

      qo.push_back(Q::zero());//I

      ostates.push_back(istates[0]);

      qo.push_back(Q(0,1));//a_down
      qo.push_back(Q(1,0));//a_up

      state.push_anni_down(i);
      ostates.push_back(state);
      state.clear();

      state.push_anni_up(i);
      ostates.push_back(state);
      state.clear();

      int row = 1;

      while(row < istates.size()){

         qo.push_back(qi[row]);//a_down
         ostates.push_back(istates[row]);

         ++row;

      }

      mpo[i].resize(Q::zero(),make_array(-qi,qp,-qp,qo));

      int column = 0;

      insert_id(mpo[i],0,0);
      insert_anni_down(mpo[i],0,1,1.0);
      insert_anni_up_s(mpo[i],0,2,1.0);

      row = 1;
      column = 3;

      //signs!
      while(row < istates.size()){

         //fermion sign!
         insert_sign(mpo[i],row,column);
         column++;
         row++;

      }

      istates = ostates;
      qi = qo;

   }

   //final occupied
   ostates.clear();
   qo.clear();

   qo.push_back(Q(0,1));//a_down
   qo.push_back(Q(1,0));//a_up

   state.push_anni_down(no - 1);
   ostates.push_back(state);
   state.clear();

   state.push_anni_up(no - 1);
   ostates.push_back(state);
   state.clear();

   int row = 1;

   while(row < istates.size()){

      qo.push_back(qi[row]);//a_down
      ostates.push_back(istates[row]);

      ++row;

   }

   mpo[no - 1].resize(Q::zero(),make_array(-qi,qp,-qp,qo));

   int column = 0;

   insert_anni_down(mpo[no - 1],0,0,1.0);
   insert_anni_up_s(mpo[no - 1],0,1,1.0);

   row = 1;
   column = 2;

   //signs!
   while(row < istates.size()){

      //fermion sign!
      insert_sign(mpo[no - 1],row,column);
      column++;
      row++;

   }

   istates = ostates;
   qi = qo;

   if(no == nv){

      ostates.clear();
      qo.clear();

      qo.push_back(Q::zero());//I
      state.push_id();
      ostates.push_back(state);
      state.clear();

      for(int i = no + 1;i < L;++i){

         qo.push_back(Q(0,1));//complementary of a^\dagger_down
         qo.push_back(Q(1,0));//complementary of a^\dagger_up

         state.push_crea_down(i);
         ostates.push_back(state);
         state.clear();

         state.push_crea_up(i);
         ostates.push_back(state);
         state.clear();

      }

      mpo[no].resize(Q::zero(),make_array(-qi,qp,-qp,qo));

      //first column is closed
      for(int row = 0;row < istates.size();++row){

         int ii = istates[row].gsite(0);
         int si = istates[row].gspin(0);

         if(si == 0)
            insert_crea_up(mpo[no],row,0,t(ii,0));
         else
            insert_crea_down_s(mpo[no],row,0,t(ii,0));

      }

      //insert the rest of the matrix: t's!
      for(int row = 0;row < istates.size();++row){

         int ii = istates[row].gsite(0);
         int si = istates[row].gspin(0);

         for(int col = 1;col < ostates.size();++col){

            int ia = ostates[col].gsite(0);
            int sa = ostates[col].gspin(0);

            if(si == sa)
               insert_sign(mpo[no],row,col,t(ii,ia - no));

         }

      }

      istates = ostates;
      qi = qo;

      //rest of the virtuals
      for(int i = no + 1;i < L - 1;++i){

         ostates.clear();
         qo.clear();

         qo.push_back(Q::zero());//I
         state.push_id();
         ostates.push_back(state);
         state.clear();


         for(int row = 3;row < istates.size();++row){

            qo.push_back(qi[row]);
            ostates.push_back(istates[row]);

         }

         mpo[i].resize(Q::zero(),make_array(-qi,qp,-qp,qo));

         //first row closed
         insert_id(mpo[i],0,0);
         insert_crea_down_s(mpo[i],1,0,1.0);
         insert_crea_up(mpo[i],2,0,1.0);

         //rest signs
         for(int row = 3;row < istates.size();++row)
            insert_sign(mpo[i],row,row-2);

         istates = ostates;
         qi = qo;

      }

      //last site
      ostates.clear();
      qo.clear();

      qo.push_back(Q::zero());//I

      mpo[L-1].resize(Q::zero(),make_array(-qi,qp,-qp,qo));

      //only closed terms
      insert_id(mpo[L-1],0,0);
      insert_crea_down_s(mpo[L-1],1,0,1.0);
      insert_crea_up(mpo[L-1],2,0,1.0);

   }
   else{//no < nv

      //first virtual: first no states are the same: so no change in ostates and qo, only id added at the end for closing terms
      state.push_id();
      ostates.push_back(state);
      state.clear();

      qo.push_back(Q::zero());//I

      mpo[no].resize(Q::zero(),make_array(-qi,qp,-qp,qo));

      //signs for the singles being copied!
      row = 0;

      while(row < istates.size()){

         //fermion sign!
         insert_sign(mpo[no],row,row);
         row++;

      }

      //insert the t's for the last column
      row = 0;
      column = istates.size();

      while(row < istates.size()){

         int j = istates[row].gsite(0);
         int sj = istates[row].gspin(0);

         //fermion sign!
         if(sj == 0)
            insert_crea_up(mpo[no],row,column,t(j,0));
         else
            insert_crea_down_s(mpo[no],row,column,t(j,0));

         ++row;

      }

      istates = ostates;
      qi = qo;

      //all virtuals until nv
      for(int i = no + 1;i < nv - 1;++i){

         int vind = i - no;

         mpo[i].resize(Q::zero(),make_array(-qi,qp,-qp,qo));

         //signs for the singles being copied!
         row = 0;

         while(row < istates.size() - 1){

            //fermion sign!
            insert_sign(mpo[i],row,row);
            row++;

         }

         //insert the t's for the last column
         row = 0;
         column = istates.size() - 1;

         while(row < istates.size() - 1){

            int j = istates[row].gsite(0);
            int sj = istates[row].gspin(0);

            //fermion sign!
            if(sj == 0)
               insert_crea_up(mpo[i],row,column,t(j,vind));
            else
               insert_crea_down_s(mpo[i],row,column,t(j,vind));

            ++row;

         }

         //last id for already closed terms
         insert_id(mpo[i],column,column);

      }

      ostates.clear();
      qo.clear();

      //switch to virtuals outgoing
      int vind = nv - no - 1;

      qo.push_back(Q::zero());//I
      state.push_id();
      ostates.push_back(state);
      state.clear();

      for(int i = nv;i < L;++i){

         qo.push_back(Q(0,1));//complementary of a^\dagger_down
         qo.push_back(Q(1,0));//complementary of a^\dagger_up

         state.push_crea_down(i);
         ostates.push_back(state);
         state.clear();

         state.push_crea_up(i);
         ostates.push_back(state);
         state.clear();

      }

      mpo[nv-1].resize(Q::zero(),make_array(-qi,qp,-qp,qo));

      //first column is closed
      for(int row = 0;row < istates.size()-1;++row){

         int ii = istates[row].gsite(0);
         int si = istates[row].gspin(0);

         if(si == 0)
            insert_crea_up(mpo[nv-1],row,0,t(ii,vind));
         else
            insert_crea_down_s(mpo[nv-1],row,0,t(ii,vind));

      }

      //insert id for the terms that are already closed
      insert_id(mpo[nv-1],istates.size() - 1,0);

      //insert the rest of the matrix: t's!
      for(int row = 0;row < istates.size();++row){

         int ii = istates[row].gsite(0);
         int si = istates[row].gspin(0);

         for(int col = 1;col < ostates.size();++col){

            int ia = ostates[col].gsite(0);
            int sa = ostates[col].gspin(0);

            if(si == sa)
               insert_sign(mpo[nv-1],row,col,t(ii,ia - no));

         }

      }

      istates = ostates;
      qi = qo;

      //rest of the virtuals
      for(int i = nv;i < L - 1;++i){

         ostates.clear();
         qo.clear();

         qo.push_back(Q::zero());//I
         state.push_id();
         ostates.push_back(state);
         state.clear();

         for(int row = 3;row < istates.size();++row){

            qo.push_back(qi[row]);
            ostates.push_back(istates[row]);

         }

         mpo[i].resize(Q::zero(),make_array(-qi,qp,-qp,qo));

         //first row closed
         insert_id(mpo[i],0,0);
         insert_crea_down_s(mpo[i],1,0,1.0);
         insert_crea_up(mpo[i],2,0,1.0);

         //rest signs
         for(int row = 3;row < istates.size();++row)
            insert_sign(mpo[i],row,row-2);

         istates = ostates;
         qi = qo;

      }

      //last site
      ostates.clear();
      qo.clear();

      qo.push_back(Q::zero());//I

      mpo[L-1].resize(Q::zero(),make_array(-qi,qp,-qp,qo));

      //only closed terms
      insert_id(mpo[L-1],0,0);
      insert_crea_down_s(mpo[L-1],1,0,1.0);
      insert_crea_up(mpo[L-1],2,0,1.0);

   }

    if(merge){

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

   }

   return mpo;

}

/**
 * @return MPO object of length L containing the T2 operator with coefficients passed through the DArray<4> object t
 */
template<class Q>
MPO<Q> T2(const DArray<4> &t,bool merge){

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

   std::vector< Ostate > ostates;
   Ostate state;

   state.push_id();
   ostates.push_back(state);
   state.clear();

   state.push_anni_down(0);
   ostates.push_back(state);
   state.clear();

   state.push_anni_up(0);
   ostates.push_back(state);
   state.clear();

   state.push_anni_up(0);
   state.push_anni_down(0);

   ostates.push_back(state);
   state.clear();

   mpo[0].resize(Q::zero(),make_array(qz,qp,-qp,qo));

   insert_id(mpo[0],0,0);
   insert_anni_down(mpo[0],0,1,1.0);
   insert_anni_up_s(mpo[0],0,2,1.0);
   insert_anni_down_anni_up(mpo[0],0,3,1.0);

   std::vector< Ostate > istates;
   istates = ostates;

   Qshapes<Quantum> qi = qo;

   if(nv < no){

      for(int i = 1;i < nv;++i){

         ostates.clear();
         qo.clear();

         qo.push_back(Q::zero());//I

         ostates.push_back(istates[0]);

         qo.push_back(Q(0,1));//a_down
         qo.push_back(Q(1,0));//a_up

         state.push_anni_down(i);
         ostates.push_back(state);
         state.clear();

         state.push_anni_up(i);
         ostates.push_back(state);
         state.clear();

         int row = 1;

         while(istates[row].size() == 1){

            qo.push_back(qi[row]);//a_down
            ostates.push_back(istates[row]);

            ++row;

         }

         qo.push_back(Q(1,1));//a_up a_down

         state.push_anni_up(i);
         state.push_anni_down(i);
         ostates.push_back(state);
         state.clear();

         //add a down
         row = 1;

         while(istates[row].size() == 1){

            state.push_anni_down(i);
            state.insert(state.end(),istates[row].begin(),istates[row].end());
            ostates.push_back(state);
            state.clear();

            Quantum tmp = qi[row];
            tmp.anni_down();
            qo.push_back(tmp);

            ++row;

         }

         //add an up
         row = 1;

         while(istates[row].size() == 1){

            state.push_anni_up(i);
            state.insert(state.end(),istates[row].begin(),istates[row].end());
            ostates.push_back(state);
            state.clear();

            Quantum tmp = qi[row];
            tmp.anni_up();
            qo.push_back(tmp);

            ++row;

         }

         //id for the pairs coming in
         while(row < istates.size()){

            qo.push_back(qi[row]);
            ostates.push_back(istates[row]);
            ++row;

         }

         mpo[i].resize(Q::zero(),make_array(-qi,qp,-qp,qo));

         int column = 0;

         insert_id(mpo[i],0,0);
         insert_anni_down(mpo[i],0,1,1.0);
         insert_anni_up_s(mpo[i],0,2,1.0);

         row = 1;
         column = 3;

         //signs!
         while(istates[row].size() == 1){

            //fermion sign!
            insert_sign(mpo[i],row,column);
            column++;
            row++;

         }

         //double remove
         insert_anni_down_anni_up(mpo[i],0,column,1.0);
         column++;

         //anni down
         row = 1;

         while(istates[row].size() == 1){

            insert_anni_down_s(mpo[i],row,column,1.0);
            column++;
            row++;

         }

         //anni up
         row = 1;

         while(istates[row].size() == 1){

            //fermion sign!
            insert_anni_up(mpo[i],row,column,1.0);
            ++column;
            ++row;

         }

         while(row < istates.size()){

            //finally identity
            insert_id(mpo[i],row,column);
            ++row;
            ++column;

         }

         istates = ostates;
         qi = qo;

      }

      //switch from incoming to outgoing doubles
      ostates.clear();
      qo.clear();

      Ostate istate;

      for(int i = no;i < L;++i){

         //first up
         istate.push_crea_up(i);

         state = istate;
         state.push_crea_down(i);
         ostates.push_back(state);
         state.clear();

         qo.push_back(Q(1,1));

         for(int j = i + 1;j < L;++j){

            //up up
            state = istate;
            state.push_crea_up(j);
            ostates.push_back(state);
            state.clear();

            qo.push_back(Q(2,0));

            //up down
            state = istate;
            state.push_crea_down(j);
            ostates.push_back(state);
            state.clear();

            qo.push_back(Q(1,1));

         }

         istate.clear();

         //first down
         istate.push_crea_down(i);

         for(int j = i + 1;j < L;++j){

            //down up
            state = istate;
            state.push_crea_up(j);
            ostates.push_back(state);
            state.clear();

            qo.push_back(Q(1,1));

            //down down
            state = istate;
            state.push_crea_down(j);
            ostates.push_back(state);
            state.clear();

            qo.push_back(Q(0,2));

         }

         istate.clear();

      }

      //then the incoming singles
      qo.push_back(Q(0,1));//a_down
      qo.push_back(Q(1,0));//a_up

      state.push_anni_down(nv);
      ostates.push_back(state);
      state.clear();

      state.push_anni_up(nv);
      ostates.push_back(state);
      state.clear();

      int row = 1;

      while(istates[row].size() == 1){

         qo.push_back(qi[row]);//a_down
         ostates.push_back(istates[row]);

         ++row;

      }

      //finally the id
      qo.push_back(Q::zero());//I
      ostates.push_back(istates[0]);

      mpo[nv].resize(Q::zero(),make_array(-qi,qp,-qp,qo));

      //first row
      int col = 0;

      while(ostates[col].size() == 2){

         int ai = ostates[col].gsite(0);
         int as = ostates[col].gspin(0);

         int bi = ostates[col].gsite(1);
         int bs = ostates[col].gspin(1);

         if(as == 0 && bs == 1)//- because of extra minus sign!
            insert_anni_down_anni_up(mpo[nv],0,col,-t(nv,nv,ai - no,bi - no));
         else if(as == 1 && bs == 0)
            insert_anni_down_anni_up(mpo[nv],0,col,t(nv,nv,ai - no,bi - no));

         ++col;

      }

      //insert single for 'outgoing' incoming singles
      while(col < ostates.size() - 1){

         int ci = ostates[col].gsite(0);
         int cs = ostates[col].gspin(0);

         if(ci == nv && cs == 0)//insert anni_up
            insert_anni_up_s(mpo[nv],0,col,1.0);
         else if(ci == nv && cs == 1)
            insert_anni_down(mpo[nv],0,col,1.0);

         ++col;

      }

      //insert id finally
      insert_id(mpo[nv],0,col);

      //singles coming in 
      row = 1;

      while(istates[row].size() == 1){

         //to outgoing doubles
         col = 0;

         while(ostates[col].size() == 2){

            double val;

            int op = Ostate::get_single_complement_T2_bis(nv,istates[row],ostates[col],t,val);

            if(op == 0)
               insert_anni_up(mpo[nv],row,col,val);
            else if(op == 1)
               insert_anni_down_s(mpo[nv],row,col,val);

            ++col;

         }

         //to incoming singles: add sign
         while(col < ostates.size() - 1){

            if(istates[row] == ostates[col])
               insert_sign(mpo[nv],row,col);

            ++col;

         }

         ++row;

      }

      //doubles coming in to doubles going out
      while(row < istates.size()){

         col = 0;

         while(ostates[col].size() == 2){

            double val;

            int op = Ostate::get_double_complement_T2(istates[row],ostates[col],t,val);

            if(op == 1)
               insert_id(mpo[nv],row,col,val);

            ++col;

         }

         ++row;

      }

      istates = ostates;
      qi = qo;

      for(int i = nv + 1;i < no - 1;++i){

         ostates.clear();
         qo.clear();

         row = 0;

         //outgoing doubles stay the same
         while(istates[row].size() == 2){

            ostates.push_back(istates[row]);
            qo.push_back(qi[row]);

            ++row;

         }

         //then the incoming singles
         qo.push_back(Q(0,1));//a_down
         qo.push_back(Q(1,0));//a_up

         state.push_anni_down(i);
         ostates.push_back(state);
         state.clear();

         state.push_anni_up(i);
         ostates.push_back(state);
         state.clear();

         while(row < istates.size() - 1){

            qo.push_back(qi[row]);//a_down
            ostates.push_back(istates[row]);

            ++row;

         }

         //end with unity
         qo.push_back(Q::zero());//I
         ostates.push_back(istates[istates.size() - 1]);

         mpo[i].resize(Q::zero(),make_array(-qi,qp,-qp,qo));

         row = 0;

         //incoming 'outgoing' doubles
         while(istates[row].size() == 2){

            col = 0;

            while(ostates[col].size() == 2){

               if(istates[row] == ostates[col])
                  insert_id(mpo[i],row,col);

               ++col;

            }

            ++row;

         }

         //incoming singles
         while(row < istates.size() - 1){

            //to outgoing doubles
            col = 0;

            while(ostates[col].size() == 2){

               double val;

               int op = Ostate::get_single_complement_T2_bis(i,istates[row],ostates[col],t,val);

               if(op == 0)
                  insert_anni_up(mpo[i],row,col,val);
               else if(op == 1)
                  insert_anni_down_s(mpo[i],row,col,val);

               ++col;

            }

            //to incoming singles: add sign
            while(col < ostates.size() - 1){

               if(istates[row] == ostates[col])
                  insert_sign(mpo[i],row,col);

               ++col;

            }

            ++row;

         }

         //finally incoming id
         int col = 0;

         while(ostates[col].size() == 2){

            int ai = ostates[col].gsite(0);
            int as = ostates[col].gspin(0);

            int bi = ostates[col].gsite(1);
            int bs = ostates[col].gspin(1);

            if(as == 0 && bs == 1)//- because of extra minus sign!
               insert_anni_down_anni_up(mpo[i],istates.size() - 1,col,-t(i,i,ai - no,bi - no));
            else if(as == 1 && bs == 0)
               insert_anni_down_anni_up(mpo[i],istates.size() - 1,col,t(i,i,ai - no,bi - no));

            ++col;

         }

         //insert single for 'outgoing' incoming singles
         while(col < ostates.size() - 1){

            int ci = ostates[col].gsite(0);
            int cs = ostates[col].gspin(0);

            if(ci == i && cs == 0)//insert anni_up
               insert_anni_up_s(mpo[i],istates.size() - 1,col,1.0);
            else if(ci == i && cs == 1)
               insert_anni_down(mpo[i],istates.size() - 1,col,1.0);

            ++col;

         }

         //insert id finally
         insert_id(mpo[i],istates.size() - 1,ostates.size() - 1);

         qi = qo;
         istates = ostates;

      }

      //site no  - 1: close down all the singles
      ostates.clear();
      qo.clear();

      row = 0;

      //only outgoing doubles
      while(istates[row].size() == 2){

         ostates.push_back(istates[row]);
         qo.push_back(qi[row]);

         ++row;

      }

      mpo[no - 1].resize(Q::zero(),make_array(-qi,qp,-qp,qo));

      row = 0;

      //incoming 'outgoing' doubles
      while(istates[row].size() == 2){

         col = 0;

         while(ostates[col].size() == 2){

            if(istates[row] == ostates[col])
               insert_id(mpo[no - 1],row,col);

            ++col;

         }

         ++row;

      }

      //incoming singles
      while(row < istates.size() - 1){

         //to outgoing doubles
         col = 0;

         while(ostates[col].size() == 2){

            double val;

            int op = Ostate::get_single_complement_T2_bis(no - 1,istates[row],ostates[col],t,val);

            if(op == 0)
               insert_anni_up(mpo[no - 1],row,col,val);
            else if(op == 1)
               insert_anni_down_s(mpo[no - 1],row,col,val);

            ++col;

         }

         ++row;

      }

      //finally incoming id
      col = 0;

      while(ostates[col].size() == 2){

         int ai = ostates[col].gsite(0);
         int as = ostates[col].gspin(0);

         int bi = ostates[col].gsite(1);
         int bs = ostates[col].gspin(1);

         if(as == 0 && bs == 1)//- because of extra minus sign!
            insert_anni_down_anni_up(mpo[no - 1],istates.size() - 1,col,-t(no - 1,no - 1,ai - no,bi - no));
         else if(as == 1 && bs == 0)
            insert_anni_down_anni_up(mpo[no - 1],istates.size() - 1,col,t(no - 1,no - 1,ai - no,bi - no));

         ++col;

      }

      istates = ostates;
      qi = qo;

      //first virtual: i = no
      qo.clear();
      ostates.clear();

      //start with id for closed states
      state.push_id();
      ostates.push_back(state);
      state.clear();

      qo.push_back(Q::zero());

      //going to single operators
      for(int j = no + 1;j < L;++j){

         state.push_crea_up(j);
         ostates.push_back(state);
         state.clear();

         qo.push_back(Q(1,0));

         state.push_crea_down(j);
         ostates.push_back(state);
         state.clear();

         qo.push_back(Q(0,1));

      }

      //going to double operators: copy from incoming
      for(int row = 0;row < istates.size();++row){

         if(istates[row].gsite(0) != no){

            ostates.push_back(istates[row]);
            qo.push_back(qi[row]);

         }

      }

      mpo[no].resize(Q::zero(),make_array(-qi,qp,-qp,qo));

      //first column close with pairs
      for(int row = 0;row < istates.size();++row){

         if(istates[row].gsite(0) == no && istates[row].gsite(1) == no)
            insert_crea_up_crea_down(mpo[no],row,0,1.0);

      }

      //singles going out, add one operator
      col = 1;

      while(ostates[col].size() == 1){

         int bi = ostates[col].gsite(0);
         int bs = ostates[col].gspin(0);

         for(int row = 0;row < istates.size();++row){

            if(bi == istates[row].gsite(1) && bs == istates[row].gspin(1)){//same outgoing state

               if(istates[row].gsite(0) == no){//first incoming state is site

                  if(istates[row].gspin(0) == 0)
                     insert_crea_up_s(mpo[no],row,col,1.0);
                  else
                     insert_crea_down(mpo[no],row,col,1.0);

               }

            }

         }

         ++col;

      }

      //doubles going out insert id
      while(col < ostates.size()){

         int cai = ostates[col].gsite(0);
         int cas = ostates[col].gspin(0);

         int cbi = ostates[col].gsite(1);
         int cbs = ostates[col].gspin(1);

         for(int row = 0;row < istates.size();++row){

            int rai = istates[row].gsite(0);
            int ras = istates[row].gspin(0);

            int rbi = istates[row].gsite(1);
            int rbs = istates[row].gspin(1);

            if(cai == rai && cas == ras && cbi == rbi && cbs == rbs)
               insert_id(mpo[no],row,col);

         }

         ++col;

      }

      istates = ostates;
      qi = qo;

      //rest of the virtuals
      for(int i = no + 1;i < L - 1;++i){

         qo.clear();
         ostates.clear();

         //start with id for closed states
         state.push_id();
         ostates.push_back(state);
         state.clear();

         qo.push_back(Q::zero());

         //going to single operators: copy from incoming
         int row = 1;

         while(istates[row].size() == 1){

            if(istates[row].gsite(0) != i){

               ostates.push_back(istates[row]);
               qo.push_back(qi[row]);

            }

            ++row;

         }

         //going to double operators: copy from incoming
         while(row < istates.size()){

            if(istates[row].gsite(0) != i){

               ostates.push_back(istates[row]);
               qo.push_back(qi[row]);

            }

            ++row;

         }

         mpo[i].resize(Q::zero(),make_array(-qi,qp,-qp,qo));

         //first column closed
         insert_id(mpo[i],0,0);

         //close the singles
         row = 1;

         while(istates[row].size() == 1){

            if(istates[row].gsite(0) == i){

               if(istates[row].gspin(0) == 0)
                  insert_crea_up(mpo[i],row,0,1.0);
               else
                  insert_crea_down_s(mpo[i],row,0,1.0);

            }

            ++row;

         }

         //close the doubles
         while(row < istates.size()){

            if(istates[row].gsite(0) == i && istates[row].gsite(1) == i)
               insert_crea_up_crea_down(mpo[i],row,0,1.0);

            ++row;

         }

         //next the single outgoing columns
         int col = 1;

         while(ostates[col].size() == 1){

            int ci = ostates[col].gsite(0);
            int cs = ostates[col].gspin(0);

            //insert sign for incoming singles
            row = 1;

            while(istates[row].size() == 1){

               int ri = istates[row].gsite(0);
               int rs = istates[row].gspin(0);

               if(ri == ci && rs == cs)
                  insert_sign(mpo[i],row,col);

               ++row;

            }

            //insert operator for incoming doubles
            while(row < istates.size()){

               if(ci == istates[row].gsite(1) && cs == istates[row].gspin(1)){//same outgoing state

                  if(istates[row].gsite(0) == i){//first incoming state is site

                     if(istates[row].gspin(0) == 0)
                        insert_crea_up_s(mpo[i],row,col,1.0);
                     else
                        insert_crea_down(mpo[i],row,col,1.0);

                  }

               }

               ++row;

            }

            ++col;

         }

         //find row where incoming doubles begin
         int rd = 1;

         while(istates[rd].size() == 1)
            ++rd;

         //insert id for doubles coming in and going out
         while(col < ostates.size()){

            int cai = ostates[col].gsite(0);
            int cas = ostates[col].gspin(0);

            int cbi = ostates[col].gsite(1);
            int cbs = ostates[col].gspin(1);

            for(int row = rd;row < istates.size();++row){

               int rai = istates[row].gsite(0);
               int ras = istates[row].gspin(0);

               int rbi = istates[row].gsite(1);
               int rbs = istates[row].gspin(1);

               if(cai == rai && cas == ras && cbi == rbi && cbs == rbs)
                  insert_id(mpo[i],row,col);

            }

            ++col;

         }

         istates = ostates;
         qi = qo;

      }

      //last site: only 4 coming in
      mpo[L - 1].resize(Q::zero(),make_array(-qi,qp,-qp,qz));

      insert_id(mpo[L - 1],0,0);
      insert_crea_up(mpo[L - 1],1,0,1.0);
      insert_crea_down_s(mpo[L - 1],2,0,1.0);
      insert_crea_up_crea_down(mpo[L - 1],3,0,1.0);

   }
   else{//nv >= no

      for(int i = 1;i < no - 1;++i){

         ostates.clear();
         qo.clear();

         qo.push_back(Q::zero());//I

         ostates.push_back(istates[0]);

         qo.push_back(Q(0,1));//a_down
         qo.push_back(Q(1,0));//a_up

         state.push_anni_down(i);
         ostates.push_back(state);
         state.clear();

         state.push_anni_up(i);
         ostates.push_back(state);
         state.clear();

         int row = 1;

         while(istates[row].size() == 1){

            qo.push_back(qi[row]);//a_down
            ostates.push_back(istates[row]);

            ++row;

         }

         qo.push_back(Q(1,1));//a_up a_down

         state.push_anni_up(i);
         state.push_anni_down(i);
         ostates.push_back(state);
         state.clear();

         //add a down
         row = 1;

         while(istates[row].size() == 1){

            state.push_anni_down(i);
            state.insert(state.end(),istates[row].begin(),istates[row].end());
            ostates.push_back(state);
            state.clear();

            Quantum tmp = qi[row];
            tmp.anni_down();
            qo.push_back(tmp);

            ++row;

         }

         //add an up
         row = 1;

         while(istates[row].size() == 1){

            state.push_anni_up(i);
            state.insert(state.end(),istates[row].begin(),istates[row].end());
            ostates.push_back(state);
            state.clear();

            Quantum tmp = qi[row];
            tmp.anni_up();
            qo.push_back(tmp);

            ++row;

         }

         //id for the pairs coming in
         while(row < istates.size()){

            qo.push_back(qi[row]);
            ostates.push_back(istates[row]);
            ++row;

         }

         mpo[i].resize(Q::zero(),make_array(-qi,qp,-qp,qo));

         int column = 0;

         insert_id(mpo[i],0,0);
         insert_anni_down(mpo[i],0,1,1.0);
         insert_anni_up_s(mpo[i],0,2,1.0);

         row = 1;
         column = 3;

         //signs!
         while(istates[row].size() == 1){

            //fermion sign!
            insert_sign(mpo[i],row,column);
            column++;
            row++;

         }

         //double remove
         insert_anni_down_anni_up(mpo[i],0,column,1.0);
         column++;

         //anni down
         row = 1;

         while(istates[row].size() == 1){

            insert_anni_down_s(mpo[i],row,column,1.0);
            column++;
            row++;

         }

         //anni up
         row = 1;

         while(istates[row].size() == 1){

            //fermion sign!
            insert_anni_up(mpo[i],row,column,1.0);
            ++column;
            ++row;

         }

         while(row < istates.size()){

            //finally identity
            insert_id(mpo[i],row,column);
            ++row;
            ++column;

         }

         istates = ostates;
         qi = qo;

      }

      if(no == nv){

         //last occupied: i = no - 1: switch from incoming to outgoing
         ostates.clear();
         qo.clear();

         Ostate istate;

         for(int i = no;i < L;++i){

            //first up
            istate.push_crea_up(i);

            state = istate;
            state.push_crea_down(i);
            ostates.push_back(state);
            state.clear();

            qo.push_back(Q(1,1));

            for(int j = i + 1;j < L;++j){

               //up up
               state = istate;
               state.push_crea_up(j);
               ostates.push_back(state);
               state.clear();

               qo.push_back(Q(2,0));

               //up down
               state = istate;
               state.push_crea_down(j);
               ostates.push_back(state);
               state.clear();

               qo.push_back(Q(1,1));

            }

            istate.clear();

            //first down
            istate.push_crea_down(i);

            for(int j = i + 1;j < L;++j){

               //down up
               state = istate;
               state.push_crea_up(j);
               ostates.push_back(state);
               state.clear();

               qo.push_back(Q(1,1));

               //down down
               state = istate;
               state.push_crea_down(j);
               ostates.push_back(state);
               state.clear();

               qo.push_back(Q(0,2));

            }

            istate.clear();

         }

         mpo[no-1].resize(Q::zero(),make_array(-qi,qp,-qp,qo));

         //first row
         for(int col = 0;col < ostates.size();++col){

            int ai = ostates[col].gsite(0);
            int as = ostates[col].gspin(0);

            int bi = ostates[col].gsite(1);
            int bs = ostates[col].gspin(1);

            if(as == 0 && bs == 1)//- because of extra minus sign!
               insert_anni_down_anni_up(mpo[no-1],0,col,-t(no-1,no-1,ai - no,bi - no));
            else if(as == 1 && bs == 0)
               insert_anni_down_anni_up(mpo[no-1],0,col,t(no-1,no-1,ai - no,bi - no));

         }

         //singles coming in
         int row = 1;

         while(istates[row].size() == 1){

            for(int col = 0;col < ostates.size();++col){

               double val;

               int op = Ostate::get_single_complement_T2_bis(no-1,istates[row],ostates[col],t,val);

               if(op == 0)
                  insert_anni_up(mpo[no - 1],row,col,val);
               else if(op == 1)
                  insert_anni_down_s(mpo[no - 1],row,col,val);

            }

            ++row;

         }

         //doubles coming in
         while(row < istates.size()){

            for(int col = 0;col < ostates.size();++col){

               double val;

               int op = Ostate::get_double_complement_T2(istates[row],ostates[col],t,val);

               if(op == 1)
                  insert_id(mpo[no - 1],row,col,val);

            }

            ++row;

         }

         istates = ostates;
         qi = qo;

         //first virtual: i = no
         qo.clear();
         ostates.clear();

         //start with id for closed states
         state.push_id();
         ostates.push_back(state);
         state.clear();

         qo.push_back(Q::zero());

         //going to single operators
         for(int j = no + 1;j < L;++j){

            state.push_crea_up(j);
            ostates.push_back(state);
            state.clear();

            qo.push_back(Q(1,0));

            state.push_crea_down(j);
            ostates.push_back(state);
            state.clear();

            qo.push_back(Q(0,1));

         }

         //going to double operators: copy from incoming
         for(int row = 0;row < istates.size();++row){

            if(istates[row].gsite(0) != no){

               ostates.push_back(istates[row]);
               qo.push_back(qi[row]);

            }

         }

         mpo[no].resize(Q::zero(),make_array(-qi,qp,-qp,qo));

         //first column close with pairs
         for(int row = 0;row < istates.size();++row){

            if(istates[row].gsite(0) == no && istates[row].gsite(1) == no)
               insert_crea_up_crea_down(mpo[no],row,0,1.0);

         }

         //singles going out, add one operator
         int col = 1;

         while(ostates[col].size() == 1){

            int bi = ostates[col].gsite(0);
            int bs = ostates[col].gspin(0);

            for(int row = 0;row < istates.size();++row){

               if(bi == istates[row].gsite(1) && bs == istates[row].gspin(1)){//same outgoing state

                  if(istates[row].gsite(0) == no){//first incoming state is site

                     if(istates[row].gspin(0) == 0)
                        insert_crea_up_s(mpo[no],row,col,1.0);
                     else
                        insert_crea_down(mpo[no],row,col,1.0);

                  }

               }

            }

            ++col;

         }

         //doubles going out insert id
         while(col < ostates.size()){

            int cai = ostates[col].gsite(0);
            int cas = ostates[col].gspin(0);

            int cbi = ostates[col].gsite(1);
            int cbs = ostates[col].gspin(1);

            for(int row = 0;row < istates.size();++row){

               int rai = istates[row].gsite(0);
               int ras = istates[row].gspin(0);

               int rbi = istates[row].gsite(1);
               int rbs = istates[row].gspin(1);

               if(cai == rai && cas == ras && cbi == rbi && cbs == rbs)
                  insert_id(mpo[no],row,col);

            }

            ++col;

         }

         istates = ostates;
         qi = qo;

         //rest of the virtuals
         for(int i = no + 1;i < L - 1;++i){

            qo.clear();
            ostates.clear();

            //start with id for closed states
            state.push_id();
            ostates.push_back(state);
            state.clear();

            qo.push_back(Q::zero());

            //going to single operators: copy from incoming
            int row = 1;

            while(istates[row].size() == 1){

               if(istates[row].gsite(0) != i){

                  ostates.push_back(istates[row]);
                  qo.push_back(qi[row]);

               }

               ++row;

            }

            //going to double operators: copy from incoming
            while(row < istates.size()){

               if(istates[row].gsite(0) != i){

                  ostates.push_back(istates[row]);
                  qo.push_back(qi[row]);

               }

               ++row;

            }

            mpo[i].resize(Q::zero(),make_array(-qi,qp,-qp,qo));

            //first column closed
            insert_id(mpo[i],0,0);

            //close the singles
            row = 1;

            while(istates[row].size() == 1){

               if(istates[row].gsite(0) == i){

                  if(istates[row].gspin(0) == 0)
                     insert_crea_up(mpo[i],row,0,1.0);
                  else
                     insert_crea_down_s(mpo[i],row,0,1.0);

               }

               ++row;

            }

            //close the doubles
            while(row < istates.size()){

               if(istates[row].gsite(0) == i && istates[row].gsite(1) == i)
                  insert_crea_up_crea_down(mpo[i],row,0,1.0);

               ++row;

            }

            //next the single outgoing columns
            int col = 1;

            while(ostates[col].size() == 1){

               int ci = ostates[col].gsite(0);
               int cs = ostates[col].gspin(0);

               //insert sign for incoming singles
               row = 1;

               while(istates[row].size() == 1){

                  int ri = istates[row].gsite(0);
                  int rs = istates[row].gspin(0);

                  if(ri == ci && rs == cs)
                     insert_sign(mpo[i],row,col);

                  ++row;

               }

               //insert operator for incoming doubles
               while(row < istates.size()){

                  if(ci == istates[row].gsite(1) && cs == istates[row].gspin(1)){//same outgoing state

                     if(istates[row].gsite(0) == i){//first incoming state is site

                        if(istates[row].gspin(0) == 0)
                           insert_crea_up_s(mpo[i],row,col,1.0);
                        else
                           insert_crea_down(mpo[i],row,col,1.0);

                     }

                  }

                  ++row;

               }

               ++col;

            }

            //find row where incoming doubles begin
            int rd = 1;

            while(istates[rd].size() == 1)
               ++rd;

            //insert id for doubles coming in and going out
            while(col < ostates.size()){

               int cai = ostates[col].gsite(0);
               int cas = ostates[col].gspin(0);

               int cbi = ostates[col].gsite(1);
               int cbs = ostates[col].gspin(1);

               for(int row = rd;row < istates.size();++row){

                  int rai = istates[row].gsite(0);
                  int ras = istates[row].gspin(0);

                  int rbi = istates[row].gsite(1);
                  int rbs = istates[row].gspin(1);

                  if(cai == rai && cas == ras && cbi == rbi && cbs == rbs)
                     insert_id(mpo[i],row,col);

               }

               ++col;

            }

            istates = ostates;
            qi = qo;

         }

         //last site: only 4 coming in
         mpo[L - 1].resize(Q::zero(),make_array(-qi,qp,-qp,qz));

         insert_id(mpo[L - 1],0,0);
         insert_crea_up(mpo[L - 1],1,0,1.0);
         insert_crea_down_s(mpo[L - 1],2,0,1.0);
         insert_crea_up_crea_down(mpo[L - 1],3,0,1.0);

      }
      else{

         //last occupied: i = no - 1
         ostates.clear();
         qo.clear();

         qo.push_back(Q(1,1));//a_up a_down

         state.push_anni_up(no-1);
         state.push_anni_down(no-1);
         ostates.push_back(state);
         state.clear();

         //add a down
         int row = 1;

         while(istates[row].size() == 1){

            state.push_anni_down(no-1);
            state.insert(state.end(),istates[row].begin(),istates[row].end());
            ostates.push_back(state);
            state.clear();

            Quantum tmp = qi[row];
            tmp.anni_down();
            qo.push_back(tmp);

            ++row;

         }

         //add an up
         row = 1;

         while(istates[row].size() == 1){

            state.push_anni_up(no-1);
            state.insert(state.end(),istates[row].begin(),istates[row].end());
            ostates.push_back(state);
            state.clear();

            Quantum tmp = qi[row];
            tmp.anni_up();
            qo.push_back(tmp);

            ++row;

         }

         //id for the pairs coming in
         while(row < istates.size()){

            qo.push_back(qi[row]);
            ostates.push_back(istates[row]);
            ++row;

         }

         //start filling
         mpo[no - 1].resize(Q::zero(),make_array(-qi,qp,-qp,qo));

         //double remove
         insert_anni_down_anni_up(mpo[no-1],0,0,1.0);

         //anni down
         row = 1;
         int column = 1;

         while(istates[row].size() == 1){

            insert_anni_down_s(mpo[no-1],row,column,1.0);
            column++;
            row++;

         }

         //anni up
         row = 1;

         while(istates[row].size() == 1){

            //fermion sign!
            insert_anni_up(mpo[no-1],row,column,1.0);
            ++column;
            ++row;

         }

         while(row < istates.size()){

            //finally identity
            insert_id(mpo[no-1],row,column);
            ++row;
            ++column;

         }

         istates = ostates;
         qi = qo;

         if(no + 1 == nv){

            ostates.clear();
            qo.clear();

            //first col closed
            state.push_id();
            ostates.push_back(state);
            state.clear();

            qo.push_back(Q::zero());

            //singles
            for(int i = no + 1;i < L;++i){

               state.push_crea_up(i);
               ostates.push_back(state);
               state.clear();

               qo.push_back(Q(1,0));

               state.push_crea_down(i);
               ostates.push_back(state);
               state.clear();

               qo.push_back(Q(0,1));

            }

            //pairs
            Ostate istate;

            for(int i = no + 1;i < L;++i){

               //first up
               istate.push_crea_up(i);

               state = istate;
               state.push_crea_down(i);
               ostates.push_back(state);
               state.clear();

               qo.push_back(Q(1,1));

               for(int j = i + 1;j < L;++j){

                  //up up
                  state = istate;
                  state.push_crea_up(j);
                  ostates.push_back(state);
                  state.clear();

                  qo.push_back(Q(2,0));

                  //up down
                  state = istate;
                  state.push_crea_down(j);
                  ostates.push_back(state);
                  state.clear();

                  qo.push_back(Q(1,1));

               }

               istate.clear();

               //first down
               istate.push_crea_down(i);

               for(int j = i + 1;j < L;++j){

                  //down up
                  state = istate;
                  state.push_crea_up(j);
                  ostates.push_back(state);
                  state.clear();

                  qo.push_back(Q(1,1));

                  //down down
                  state = istate;
                  state.push_crea_down(j);
                  ostates.push_back(state);
                  state.clear();

                  qo.push_back(Q(0,2));

               }

               istate.clear();

            }

            mpo[no].resize(Q::zero(),make_array(-qi,qp,-qp,qo));

            //first column closed: insert pairs
            for(int row = 0;row < istates.size();++row){

               //lets call in i,j
               int i = istates[row].gsite(1);
               int j = istates[row].gsite(0);

               int si = istates[row].gspin(1);
               int sj = istates[row].gspin(0);

               if(si == 0 && sj == 1)
                  insert_crea_up_crea_down(mpo[no],row,0,t(i,j,0,0));
               else if(si == 1 && sj == 0)
                  insert_crea_up_crea_down(mpo[no],row,0,-t(i,j,0,0));

            }

            //singles going out:
            int col = 1;

            while(ostates[col].size() == 1){

               for(int row = 0;row < istates.size();++row){

                  double val;

                  int op = Ostate::get_single_complement_T2(no,istates[row],ostates[col],t,val);

                  if(op == 0)
                     insert_crea_up_s(mpo[no],row,col,val);
                  else if(op == 1)
                     insert_crea_down(mpo[no],row,col,val);

               }

               ++col;

            }

            //switch to outgoing pairs
            while(col < ostates.size()){

               for(int row = 0;row < istates.size();++row){

                  double val;

                  int op = Ostate::get_double_complement_T2(istates[row],ostates[col],t,val);

                  if(op == 1)
                     insert_id(mpo[no],row,col,val);

               }

               ++col;

            }

            istates = ostates;
            qi = qo;

            //rest of the virtuals
            for(int i = no + 1;i < L - 1;++i){

               qo.clear();
               ostates.clear();

               //start with id for closed states
               state.push_id();
               ostates.push_back(state);
               state.clear();

               qo.push_back(Q::zero());

               //going to single operators: copy from incoming
               int row = 1;

               while(istates[row].size() == 1){

                  if(istates[row].gsite(0) != i){

                     ostates.push_back(istates[row]);
                     qo.push_back(qi[row]);

                  }

                  ++row;

               }

               //going to double operators: copy from incoming
               while(row < istates.size()){

                  if(istates[row].gsite(0) != i){

                     ostates.push_back(istates[row]);
                     qo.push_back(qi[row]);

                  }

                  ++row;

               }

               mpo[i].resize(Q::zero(),make_array(-qi,qp,-qp,qo));

               //first column closed
               insert_id(mpo[i],0,0);

               //close the singles
               row = 1;

               while(istates[row].size() == 1){

                  if(istates[row].gsite(0) == i){

                     if(istates[row].gspin(0) == 0)
                        insert_crea_up(mpo[i],row,0,1.0);
                     else
                        insert_crea_down_s(mpo[i],row,0,1.0);

                  }

                  ++row;

               }

               //close the doubles
               while(row < istates.size()){

                  if(istates[row].gsite(0) == i && istates[row].gsite(1) == i)
                     insert_crea_up_crea_down(mpo[i],row,0,1.0);

                  ++row;

               }

               //next the single outgoing columns
               int col = 1;

               while(ostates[col].size() == 1){

                  int ci = ostates[col].gsite(0);
                  int cs = ostates[col].gspin(0);

                  //insert sign for incoming singles
                  row = 1;

                  while(istates[row].size() == 1){

                     int ri = istates[row].gsite(0);
                     int rs = istates[row].gspin(0);

                     if(ri == ci && rs == cs)
                        insert_sign(mpo[i],row,col);

                     ++row;

                  }

                  //insert operator for incoming doubles
                  while(row < istates.size()){

                     if(ci == istates[row].gsite(1) && cs == istates[row].gspin(1)){//same outgoing state

                        if(istates[row].gsite(0) == i){//first incoming state is site

                           if(istates[row].gspin(0) == 0)
                              insert_crea_up_s(mpo[i],row,col,1.0);
                           else
                              insert_crea_down(mpo[i],row,col,1.0);

                        }

                     }

                     ++row;

                  }

                  ++col;

               }

               //find row where incoming doubles begin
               int rd = 1;

               while(istates[rd].size() == 1)
                  ++rd;

               //insert id for doubles coming in and going out
               while(col < ostates.size()){

                  int cai = ostates[col].gsite(0);
                  int cas = ostates[col].gspin(0);

                  int cbi = ostates[col].gsite(1);
                  int cbs = ostates[col].gspin(1);

                  for(int row = rd;row < istates.size();++row){

                     int rai = istates[row].gsite(0);
                     int ras = istates[row].gspin(0);

                     int rbi = istates[row].gsite(1);
                     int rbs = istates[row].gspin(1);

                     if(cai == rai && cas == ras && cbi == rbi && cbs == rbs)
                        insert_id(mpo[i],row,col);

                  }

                  ++col;

               }

               istates = ostates;
               qi = qo;

            }

            //last site: only 4 coming in
            mpo[L - 1].resize(Q::zero(),make_array(-qi,qp,-qp,qz));

            insert_id(mpo[L - 1],0,0);
            insert_crea_up(mpo[L - 1],1,0,1.0);
            insert_crea_down_s(mpo[L - 1],2,0,1.0);
            insert_crea_up_crea_down(mpo[L - 1],3,0,1.0);

         }
         else{//nv > no + 1

            //for once the ostates and qo from previous site are correct for start of next site

            //all the virtuals: the have the signature of the operator they are going to: but opposite q-number
            for(int j = no + 1;j < L;++j){

               qo.push_back(Q(1,0));//go to crea up

               state.push_crea_up(j);
               ostates.push_back(state);
               state.clear();

               qo.push_back(Q(0,1));//go to crea down 

               state.push_crea_down(j);
               ostates.push_back(state);
               state.clear();

            }

            //last column complete closed:
            qo.push_back(Q(0,0));
            state.push_id();
            ostates.push_back(state);
            state.clear();

            //start the filling
            mpo[no].resize(Q::zero(),make_array(-qi,qp,-qp,qo));

            //copy the pairs
            for(int i = 0;i < istates.size();++i)
               insert_id(mpo[no],i,i);

            for(int column = istates.size();column < ostates.size() - 1;++column){

               for(int row = 0;row < istates.size();++row){

                  double val;

                  int op = Ostate::get_single_complement_T2(no,istates[row],ostates[column],t,val);

                  if(op == 0)
                     insert_crea_up_s(mpo[no],row,column,val);
                  else if(op == 1)
                     insert_crea_down(mpo[no],row,column,val);

               }

            }

            //last column: insert pairs
            for(int row = 0;row < istates.size();++row){

               //lets call in i,j
               int i = istates[row].gsite(1);
               int j = istates[row].gsite(0);

               int si = istates[row].gspin(1);
               int sj = istates[row].gspin(0);

               if(si == 0 && sj == 1)
                  insert_crea_up_crea_down(mpo[no],row,ostates.size() - 1,t(i,j,0,0));
               else if(si == 1 && sj == 0)
                  insert_crea_up_crea_down(mpo[no],row,ostates.size() - 1,-t(i,j,0,0));

            }

            istates = ostates;
            qi = qo;

            //next everything until nv
            for(int i = no + 1;i < nv;++i){

               //current virtual index
               int vind = i - no;

               ostates.clear();
               qo.clear();

               //identity for the pairs
               int row = 0;

               while(istates[row].size() == 2){

                  ostates.push_back(istates[row]);
                  qo.push_back(qi[row]);
                  ++row;

               }

               //the remaining virtuals: they have the signature of the operator they are going to: but opposite q-number
               for(int j = i + 1;j < L;++j){

                  qo.push_back(Q(1,0));//go to crea up

                  state.push_crea_up(j);
                  ostates.push_back(state);
                  state.clear();

                  qo.push_back(Q(0,1));//go to crea down 

                  state.push_crea_down(j);
                  ostates.push_back(state);
                  state.clear();

               }

               //last column closed
               state.push_id();
               ostates.push_back(state);
               state.clear();
               qo.push_back(Q(0,0));

               //start filling
               mpo[i].resize(Q::zero(),make_array(-qi,qp,-qp,qo));

               row = 0;
               int column = 0;

               //fill pairs id
               while(istates[row].size() == 2){

                  insert_id(mpo[i],row,column);
                  ++row;++column;

               }

               //remaining virtuals
               while(column < ostates.size() - 1){

                  row = 0;

                  while(istates[row].size() == 2){

                     double val;

                     int op = Ostate::get_single_complement_T2(i,istates[row],ostates[column],t,val);

                     if(op == 0)
                        insert_crea_up_s(mpo[i],row,column,val);
                     else if(op == 1)
                        insert_crea_down(mpo[i],row,column,val);

                     ++row;

                  }

                  while(row < istates.size() - 1){

                     if(istates[row].gsite(0) == ostates[column].gsite(0) && istates[row].gspin(0) == ostates[column].gspin(0))
                        insert_sign(mpo[i],row,column);

                     ++row;

                  }

                  ++column;

               }

               //last column: close down
               row = 0;

               while(istates[row].size() == 2){

                  //lets call in i,j
                  int ii = istates[row].gsite(1);
                  int ij = istates[row].gsite(0);

                  int si = istates[row].gspin(1);
                  int sj = istates[row].gspin(0);

                  if(si == 0 && sj == 1)
                     insert_crea_up_crea_down(mpo[i],row,ostates.size() - 1,t(ii,ij,vind,vind));
                  else if(si == 1 && sj == 0)
                     insert_crea_up_crea_down(mpo[i],row,ostates.size() - 1,-t(ii,ij,vind,vind));

                  ++row;

               }

               //close down the complementaries:
               insert_crea_up(mpo[i],row,ostates.size() - 1,1.0);++row;
               insert_crea_down_s(mpo[i],row,ostates.size() - 1,1.0);

               //unit for the already closed parts
               insert_id(mpo[i],istates.size() - 1,ostates.size() - 1);

               istates = ostates;
               qi = qo;

            }

            ostates.clear();
            qo.clear();

            //first col closed
            state.push_id();
            ostates.push_back(state);
            state.clear();

            qo.push_back(Q::zero());

            //singles
            for(int i = nv + 1;i < L;++i){

               state.push_crea_up(i);
               ostates.push_back(state);
               state.clear();

               qo.push_back(Q(1,0));

               state.push_crea_down(i);
               ostates.push_back(state);
               state.clear();

               qo.push_back(Q(0,1));

            }

            //pairs
            Ostate istate;

            for(int i = nv + 1;i < L;++i){

               //first up
               istate.push_crea_up(i);

               state = istate;
               state.push_crea_down(i);
               ostates.push_back(state);
               state.clear();

               qo.push_back(Q(1,1));

               for(int j = i + 1;j < L;++j){

                  //up up
                  state = istate;
                  state.push_crea_up(j);
                  ostates.push_back(state);
                  state.clear();

                  qo.push_back(Q(2,0));

                  //up down
                  state = istate;
                  state.push_crea_down(j);
                  ostates.push_back(state);
                  state.clear();

                  qo.push_back(Q(1,1));

               }

               istate.clear();

               //first down
               istate.push_crea_down(i);

               for(int j = i + 1;j < L;++j){

                  //down up
                  state = istate;
                  state.push_crea_up(j);
                  ostates.push_back(state);
                  state.clear();

                  qo.push_back(Q(1,1));

                  //down down
                  state = istate;
                  state.push_crea_down(j);
                  ostates.push_back(state);
                  state.clear();

                  qo.push_back(Q(0,2));

               }

               istate.clear();

            }

            mpo[nv].resize(Q::zero(),make_array(-qi,qp,-qp,qo));

            //incoming pairs
            int row = 0;

            while(istates[row].size() == 2){

               //lets call in i,j
               int i = istates[row].gsite(1);
               int j = istates[row].gsite(0);

               int si = istates[row].gspin(1);
               int sj = istates[row].gspin(0);

               if(si == 0 && sj == 1)
                  insert_crea_up_crea_down(mpo[nv],row,0,t(i,j,nv-no,nv-no));
               else if(si == 1 && sj == 0)
                  insert_crea_up_crea_down(mpo[nv],row,0,-t(i,j,nv-no,nv-no));

               ++row;

            }

            //incoming singles
            while(row < istates.size() - 1){

               if(istates[row].gsite(0) == nv){

                  if(istates[row].gspin(0) == 0)
                     insert_crea_up(mpo[nv],row,0,1.0);
                  else
                     insert_crea_down_s(mpo[nv],row,0,1.0);

               }

               ++row;

            }

            //for closed terms coming in
            insert_id(mpo[nv],row,0);

            //next the single columns
            int col = 1;

            while(ostates[col].size() == 1){

               row = 0;

               //doubles coming in
               while(istates[row].size() == 2){

                  double val;

                  int op = Ostate::get_single_complement_T2(nv,istates[row],ostates[col],t,val);

                  if(op == 0)
                     insert_crea_up_s(mpo[nv],row,col,val);
                  else if(op == 1)
                     insert_crea_down(mpo[nv],row,col,val);

                  ++row;

               }

               //sign for incoming to outgoing singles
               while(row < istates.size() - 1){

                  if(istates[row].gsite(0) == ostates[col].gsite(0) && istates[row].gspin(0) == ostates[col].gspin(0))
                     insert_sign(mpo[nv],row,col);

                  ++row;

               }

               ++col;

            }

            while(col < ostates.size()){

               //transform from incoming to outgoing pairs
               row = 0;

               while(istates[row].size() == 2){

                  double val;

                  int op = Ostate::get_double_complement_T2(istates[row],ostates[col],t,val);

                  if(op == 1)
                     insert_id(mpo[nv],row,col,val);

                  ++row;

               }

               ++col;

            }

            istates = ostates;
            qi = qo;

            //rest of the virtuals
            for(int i = nv + 1;i < L - 1;++i){

               qo.clear();
               ostates.clear();

               //start with id for closed states
               state.push_id();
               ostates.push_back(state);
               state.clear();

               qo.push_back(Q::zero());

               //going to single operators: copy from incoming
               int row = 1;

               while(istates[row].size() == 1){

                  if(istates[row].gsite(0) != i){

                     ostates.push_back(istates[row]);
                     qo.push_back(qi[row]);

                  }

                  ++row;

               }

               //going to double operators: copy from incoming
               while(row < istates.size()){

                  if(istates[row].gsite(0) != i){

                     ostates.push_back(istates[row]);
                     qo.push_back(qi[row]);

                  }

                  ++row;

               }

               mpo[i].resize(Q::zero(),make_array(-qi,qp,-qp,qo));

               //first column closed
               insert_id(mpo[i],0,0);

               //close the singles
               row = 1;

               while(istates[row].size() == 1){

                  if(istates[row].gsite(0) == i){

                     if(istates[row].gspin(0) == 0)
                        insert_crea_up(mpo[i],row,0,1.0);
                     else
                        insert_crea_down_s(mpo[i],row,0,1.0);

                  }

                  ++row;

               }

               //close the doubles
               while(row < istates.size()){

                  if(istates[row].gsite(0) == i && istates[row].gsite(1) == i)
                     insert_crea_up_crea_down(mpo[i],row,0,1.0);

                  ++row;

               }

               //next the single outgoing columns
               int col = 1;

               while(ostates[col].size() == 1){

                  int ci = ostates[col].gsite(0);
                  int cs = ostates[col].gspin(0);

                  //insert sign for incoming singles
                  row = 1;

                  while(istates[row].size() == 1){

                     int ri = istates[row].gsite(0);
                     int rs = istates[row].gspin(0);

                     if(ri == ci && rs == cs)
                        insert_sign(mpo[i],row,col);

                     ++row;

                  }

                  //insert operator for incoming doubles
                  while(row < istates.size()){

                     if(ci == istates[row].gsite(1) && cs == istates[row].gspin(1)){//same outgoing state

                        if(istates[row].gsite(0) == i){//first incoming state is site

                           if(istates[row].gspin(0) == 0)
                              insert_crea_up_s(mpo[i],row,col,1.0);
                           else
                              insert_crea_down(mpo[i],row,col,1.0);

                        }

                     }

                     ++row;

                  }

                  ++col;

               }

               //find row where incoming doubles begin
               int rd = 1;

               while(istates[rd].size() == 1)
                  ++rd;

               //insert id for doubles coming in and going out
               while(col < ostates.size()){

                  int cai = ostates[col].gsite(0);
                  int cas = ostates[col].gspin(0);

                  int cbi = ostates[col].gsite(1);
                  int cbs = ostates[col].gspin(1);

                  for(int row = rd;row < istates.size();++row){

                     int rai = istates[row].gsite(0);
                     int ras = istates[row].gspin(0);

                     int rbi = istates[row].gsite(1);
                     int rbs = istates[row].gspin(1);

                     if(cai == rai && cas == ras && cbi == rbi && cbs == rbs)
                        insert_id(mpo[i],row,col);

                  }

                  ++col;

               }

               istates = ostates;
               qi = qo;

            }

            //last site: only 4 coming in
            mpo[L - 1].resize(Q::zero(),make_array(-qi,qp,-qp,qz));

            insert_id(mpo[L - 1],0,0);
            insert_crea_up(mpo[L - 1],1,0,1.0);
            insert_crea_down_s(mpo[L - 1],2,0,1.0);
            insert_crea_up_crea_down(mpo[L - 1],3,0,1.0);

         }

      }

   }

   if(merge){

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

   }

   return mpo;

}

/**
 * @return MPO object of length L containing a general one-body operator \sum_ab t_ab a^+_as a_bs
 */
template<class Q>
MPO<Q> one_body(const DArray<2> &t,bool merge){

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

   mpo[0].resize(Q::zero(),make_array(qz,qp,-qp,qo));

   std::vector< Ostate > ostates;

   //identity
   Ostate state;
   state.push_id();
   ostates.push_back(state);
   state.clear();

   //a^+_up (-1)^n_down
   state.push_crea_up(0);
   ostates.push_back(state);
   state.clear();

   //a^dagger_down 
   state.push_crea_down(0);
   ostates.push_back(state);
   state.clear();

   //a_up (-1)^n_down
   state.push_anni_up(0);
   ostates.push_back(state);
   state.clear();

   //a_down
   state.push_anni_down(0);
   ostates.push_back(state);
   state.clear();

   //local term
   state.push_id();
   ostates.push_back(state);
   state.clear();

   insert_id(mpo[0],0,0);
   insert_crea_up_s(mpo[0],0,1,1.0);
   insert_crea_down(mpo[0],0,2,1.0);
   insert_anni_up_s(mpo[0],0,3,1.0);
   insert_anni_down(mpo[0],0,4,1.0);
   insert_local_ob(mpo[0],0,5,t(0,0));

   std::vector< Ostate > istates;
   istates = ostates;

   Qshapes<Quantum> qi;
   qi = qo;

   //all middle tensors
   for(int i = 1;i < L/2;++i){

      ostates.clear();
      qo.clear();

      qo.push_back(Q::zero());//I

      //identity
      state.push_id();
      ostates.push_back(state);
      state.clear();

      qo.push_back(Q(-1,0));//a^+_up

      //first create spin up
      state.push_crea_up(i);
      ostates.push_back(state);
      state.clear();

      qo.push_back(Q(0,-1));//a^+_down

      state.push_crea_down(i);
      ostates.push_back(state);
      state.clear();

      qo.push_back(Q(1,0));//a_up

      state.push_anni_up(i);
      ostates.push_back(state);
      state.clear();

      qo.push_back(Q(0,1));//a_down

      state.push_anni_down(i);
      ostates.push_back(state);
      state.clear();

      for(int j = 1;j < istates.size() - 1;++j){

         qo.push_back(qi[j]);
         ostates.push_back(istates[j]);

      }

      //last column is closed
      qo.push_back(Q::zero());//a_down
      state.push_id();
      ostates.push_back(state);
      state.clear();

      //fill the mpo!
      mpo[i].resize(Q::zero(),make_array(-qi,qp,-qp,qo));

      int row = 0;
      int column = 0;

      //identity
      insert_id(mpo[i],0,0);
      insert_crea_up_s(mpo[i],0,1,1.0);
      insert_crea_down(mpo[i],0,2,1.0);
      insert_anni_up_s(mpo[i],0,3,1.0);
      insert_anni_down(mpo[i],0,4,1.0);

      row = 1;
      column = 5;

      //signs
      for(int j = 1;j < istates.size() - 1;++j){

         insert_sign(mpo[i],row,column);
         ++row;
         ++column;

      }

      //last column: first local term
      insert_local_ob(mpo[i],0,column,t(i,i));

      row = 1;

      while(row < istates.size() - 1){

         std::vector<int> v;
         double val;

         v = Ostate::get_closing_single(i,istates[row],t,val);

         if(v[0] == 0 && v[1] == 0)//crea up
            insert_crea_up(mpo[i],row,column,val);
         else if(v[0] == 1 && v[1] == 0)//crea down
            insert_crea_down_s(mpo[i],row,column,val);
         else if(v[0] == 0 && v[1] == 1)//anni up
            insert_anni_up(mpo[i],row,column,val);
         else
            insert_anni_down_s(mpo[i],row,column,val);

         ++row;

      }

      //finally identity
      insert_id(mpo[i],row,column);

      istates = ostates;
      qi = qo;

   }

   ostates.clear();
   qo.clear();

   qo.push_back(Q::zero());//I

   //identity
   state.push_id();
   ostates.push_back(state);
   state.clear();

   //outgoing states become states where we are going to:
   for(int j = L/2 + 1;j < L;++j){

      //go to anni up
      state.push_anni_up(j);
      ostates.push_back(state);
      state.clear();

      qo.push_back(Q(-1,0));//complementary of a_up

      state.push_anni_down(j);
      ostates.push_back(state);
      state.clear();

      qo.push_back(Q(0,-1));//complementary a_down

      state.push_crea_up(j);
      ostates.push_back(state);
      state.clear();

      qo.push_back(Q(1,0));//complementary of a^+_up

      state.push_crea_down(j);
      ostates.push_back(state);
      state.clear();

      qo.push_back(Q(0,1));//complementary of a^+_down

   }

   //finally close with an id column
   state.push_id();
   ostates.push_back(state);
   state.clear();

   qo.push_back(Q::zero());

   //fill the mpo!
   mpo[L/2].resize(Q::zero(),make_array(-qi,qp,-qp,qo));

   //first column is closed
   insert_local_ob(mpo[L/2],0,0,t(L/2,L/2));

   int row = 1;

   while(row < istates.size() - 1){

      std::vector<int> v;
      double val;

      v = Ostate::get_closing_single(L/2,istates[row],t,val);

      if(v[0] == 0 && v[1] == 0)//crea up
         insert_crea_up(mpo[L/2],row,0,val);
      else if(v[0] == 1 && v[1] == 0)//crea down
         insert_crea_down_s(mpo[L/2],row,0,val);
      else if(v[0] == 0 && v[1] == 1)//anni up
         insert_anni_up(mpo[L/2],row,0,val);
      else
         insert_anni_down_s(mpo[L/2],row,0,val);

      ++row;

   }

   //finally identity for already closed terms
   insert_id(mpo[L/2],row,0);

   //first row inserts creation and annihilation operators
   row = 0;

   for(int col = 1;col < ostates.size() - 1;++col){

      int ci = ostates[col].gsite(0);
      int cs = ostates[col].gspin(0);
      int ca = ostates[col].gact(0);

      if(cs == 0 && ca == 1)
         insert_crea_up_s(mpo[L/2],0,col,t(L/2,ci));
      else if(cs == 1 && ca == 1)
         insert_crea_down(mpo[L/2],0,col,t(L/2,ci));
      else if(cs == 0 && ca == 0)
         insert_anni_up_s(mpo[L/2],0,col,t(ci,L/2));
      else
         insert_anni_down(mpo[L/2],0,col,t(ci,L/2));

   }

   //rest transfers incoming to outgoing
   for(int row = 1;row < istates.size() - 1;++row){

      int ri = istates[row].gsite(0);
      int rs = istates[row].gspin(0);
      int ra = istates[row].gact(0);

      for(int col = 1;col < ostates.size() - 1;++col){

         int ci = ostates[col].gsite(0);
         int cs = ostates[col].gspin(0);
         int ca = ostates[col].gact(0);

         if(rs == cs && ra != ca)
            insert_sign(mpo[L/2],row,col,t(ri,ci));

      }
   }

   //id for last column
   insert_id(mpo[L/2],0,ostates.size() - 1);

   istates = ostates;
   qi = qo;

   //second half goes down
   for(int i = L/2 + 1;i < L - 1;++i){

      ostates.clear();
      qo.clear();

      qo.push_back(Q::zero());//I

      //identity
      state.push_id();
      ostates.push_back(state);
      state.clear();

      //outgoing states
      for(int j = i + 1;j < L;++j){

         //go to anni up
         state.push_anni_up(j);
         ostates.push_back(state);
         state.clear();

         qo.push_back(Q(-1,0));//complementary of a_up

         state.push_anni_down(j);
         ostates.push_back(state);
         state.clear();

         qo.push_back(Q(0,-1));//complementary a_down

         state.push_crea_up(j);
         ostates.push_back(state);
         state.clear();

         qo.push_back(Q(1,0));//complementary of a^+_up

         state.push_crea_down(j);
         ostates.push_back(state);
         state.clear();

         qo.push_back(Q(0,1));//complementary of a^+_down

      }

      //finally close with an id column
      state.push_id();
      ostates.push_back(state);
      state.clear();

      qo.push_back(Q::zero());

      //fill the mpo!
      mpo[i].resize(Q::zero(),make_array(-qi,qp,-qp,qo));

      //identity for already closed terms
      insert_id(mpo[i],0,0);

      //insert correct operators for incoming complementaries
      for(int row = 1;row < istates.size() - 1;++row){

         int ri = istates[row].gsite(0);
         int rs = istates[row].gspin(0);
         int ra = istates[row].gact(0);

         if(ri == i){

            if(rs == 0 && ra == 0)//crea up
               insert_crea_up(mpo[i],row,0,1.0);
            else if(rs == 1 && ra == 0)//crea down
               insert_crea_down_s(mpo[i],row,0,1.0);
            else if(rs == 0 && ra == 1)//anni up
               insert_anni_up(mpo[i],row,0,1.0);
            else
               insert_anni_down_s(mpo[i],row,0,1.0);

         }

      }

      //insert local term
      insert_local_ob(mpo[i],istates.size() - 1,0,t(i,i));

      //signs!
      for(int row = 1;row < istates.size() - 1;++row){

         int ri = istates[row].gsite(0);
         int rs = istates[row].gspin(0);
         int ra = istates[row].gact(0);

         for(int col = 1;col < ostates.size() - 1;++col){

            int ci = ostates[col].gsite(0);
            int cs = ostates[col].gspin(0);
            int ca = ostates[col].gact(0);

            if(ri == ci && rs == cs && ra == ca)
               insert_sign(mpo[i],row,col);

         }
      }

      //last row, creation and annihilation operators
      int row = istates.size() - 1;

      for(int col = 1;col < ostates.size() - 1;++col){

         int ci = ostates[col].gsite(0);
         int cs = ostates[col].gspin(0);
         int ca = ostates[col].gact(0);

         if(cs == 0 && ca == 1)
            insert_crea_up_s(mpo[i],row,col,t(i,ci));
         else if(cs == 1 && ca == 1)
            insert_crea_down(mpo[i],row,col,t(i,ci));
         else if(cs == 0 && ca == 0)
            insert_anni_up_s(mpo[i],row,col,t(ci,i));
         else
            insert_anni_down(mpo[i],row,col,t(ci,i));

      }

      //finally insert id
      insert_id(mpo[i],istates.size() - 1,ostates.size() - 1);

      istates = ostates;
      qi = qo;

   }

   ostates.clear();
   qo.clear();

   //last site:
   qo.push_back(Q::zero());
   state.push_id();
   ostates.push_back(state);
   state.clear();

   mpo[L-1].resize(Q::zero(),make_array(-qi,qp,-qp,qo));

   //identity for already closed terms
   insert_id(mpo[L-1],0,0);

   //insert correct operators for incoming complementaries
   for(int row = 1;row < istates.size() - 1;++row){

      int rs = istates[row].gspin(0);
      int ra = istates[row].gact(0);

      if(rs == 0 && ra == 0)//crea up
         insert_crea_up(mpo[L-1],row,0,1.0);
      else if(rs == 1 && ra == 0)//crea down
         insert_crea_down_s(mpo[L-1],row,0,1.0);
      else if(rs == 0 && ra == 1)//anni up
         insert_anni_up(mpo[L-1],row,0,1.0);
      else
         insert_anni_down_s(mpo[L-1],row,0,1.0);

   }

   //insert local term
   insert_local_ob(mpo[L-1],istates.size() - 1,0,t(L-1,L-1));

   if(merge){

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

   }

   return mpo;

}

/**
 * insert identity operator in mpo O
 */
void insert_id(QSDArray<4> &O,int row,int column,double value){

   DArray<4> Ip(1,1,1,1);
   Ip = value;

   O.insert(shape(row,0,0,column),Ip);
   O.insert(shape(row,1,1,column),Ip);
   O.insert(shape(row,2,2,column),Ip);
   O.insert(shape(row,3,3,column),Ip);

}

/**
 * insert zero
 */
void insert_zero(QSDArray<4> &O,int row,int column){

   DArray<4> Ip(1,1,1,1);
   Ip = 0;

   O.insert(shape(row,0,0,column),Ip);
   O.insert(shape(row,1,1,column),Ip);
   O.insert(shape(row,2,2,column),Ip);
   O.insert(shape(row,3,3,column),Ip);

}


/**
 * insert identity operator with fermion sign
 */
void insert_sign(QSDArray<4> &O,int row,int column,double value){

   DArray<4> Ip(1,1,1,1);

   Ip = value;
   O.insert(shape(row,0,0,column),Ip);
   O.insert(shape(row,3,3,column),Ip);

   Ip = -value;
   O.insert(shape(row,1,1,column),Ip);
   O.insert(shape(row,2,2,column),Ip);

}

/**
 * insert creator of up spin
 */
void insert_crea_up(QSDArray<4> &O,int row,int column,double value){

   DArray<4> Ip(1,1,1,1);
   Ip = value;

   O.insert(shape(row,1,0,column),Ip);
   O.insert(shape(row,3,2,column),Ip);

}

/**
 * insert creator of up spin with sign down
 */
void insert_crea_up_s(QSDArray<4> &O,int row,int column,double value){

   DArray<4> Ip(1,1,1,1);
   Ip = value;

   //a^+_up (-1)^n_down
   O.insert(shape(row,1,0,column),Ip);

   Ip = -value;

   O.insert(shape(row,3,2,column),Ip);

}

/**
 * insert creator of down spin
 */
void insert_crea_down(QSDArray<4> &O,int row,int column,double value){

   DArray<4> Ip(1,1,1,1);
   Ip = value;

   //a^dagger_down 
   O.insert(shape(row,2,0,column),Ip);
   O.insert(shape(row,3,1,column),Ip);

}

/**
 * insert creator of down spin with up spin sign
 */
void insert_crea_down_s(QSDArray<4> &O,int row,int column,double value){

   DArray<4> Ip(1,1,1,1);

   //a^dagger_down 
   Ip = value;
   O.insert(shape(row,2,0,column),Ip);

   Ip = -value;
   O.insert(shape(row,3,1,column),Ip);

}

/**
 * insert annihilator of up spin
 */
void insert_anni_up(QSDArray<4> &O,int row,int column,double value){

   DArray<4> Ip(1,1,1,1);
   Ip = value;

   //a_up (-1)^n_down
   O.insert(shape(row,0,1,column),Ip);
   O.insert(shape(row,2,3,column),Ip);

}

/**
 * insert annihilator of up spin with down sign
 */
void insert_anni_up_s(QSDArray<4> &O,int row,int column,double value){

   DArray<4> Ip(1,1,1,1);
   Ip = value;

   //a_up (-1)^n_down
   O.insert(shape(row,0,1,column),Ip);

   Ip = -value;

   O.insert(shape(row,2,3,column),Ip);

}

/**
 * insert annihilator of down spin
 */
void insert_anni_down(QSDArray<4> &O,int row,int column,double value){

   DArray<4> Ip(1,1,1,1);
   Ip = value;

   //a_down
   O.insert(shape(row,0,2,column),Ip);
   O.insert(shape(row,1,3,column),Ip);

}

/**
 * insert annihilator of down spin with sign for up spin
 */
void insert_anni_down_s(QSDArray<4> &O,int row,int column,double value){

   DArray<4> Ip(1,1,1,1);

   //a_down
   Ip = value;
   O.insert(shape(row,0,2,column),Ip);

   Ip = -value;
   O.insert(shape(row,1,3,column),Ip);

}

/**
 * insert creator of an up down pair on site
 */
void insert_crea_up_crea_down(QSDArray<4> &O,int row,int column,double value){

   DArray<4> Ip(1,1,1,1);
   Ip = value;

   //a_down
   O.insert(shape(row,3,0,column),Ip);

}

/**
 * insert n_up operator
 */
void insert_crea_up_anni_up(QSDArray<4> &O,int row,int column,double value){

   DArray<4> Ip(1,1,1,1);
   Ip = value;

   O.insert(shape(row,1,1,column),Ip);
   O.insert(shape(row,3,3,column),Ip);

}

/**
 * insert create up annihilate down
 */
void insert_crea_up_anni_down(QSDArray<4> &O,int row,int column,double value){

   DArray<4> Ip(1,1,1,1);
   Ip = value;

   O.insert(shape(row,1,2,column),Ip);

}

/**
 * insert create down annihilate up
 */
void insert_crea_down_anni_up(QSDArray<4> &O,int row,int column,double value){

   DArray<4> Ip(1,1,1,1);
   Ip = value;

   O.insert(shape(row,2,1,column),Ip);

}

/**
 * insert create down annihilate down --> n_down
 */
void insert_crea_down_anni_down(QSDArray<4> &O,int row,int column,double value){

   DArray<4> Ip(1,1,1,1);
   Ip = value;

   O.insert(shape(row,2,2,column),Ip);
   O.insert(shape(row,3,3,column),Ip);

}

/**
 * insert annihilate up annihilate down: extra minus sign!
 */
void insert_anni_down_anni_up(QSDArray<4> &O,int row,int column,double value){

   DArray<4> Ip(1,1,1,1);
   Ip = -value;

   O.insert(shape(row,0,3,column),Ip);

}

/**
 * insert complementary operator for a_up: behaves as a creator of up
 */
void insert_triple_anni_up_first(QSDArray<4> &O,int row,int column,double tval,double Vval){

   DArray<4> Ip(1,1,1,1);

   Ip = tval;
   O.insert(shape(row,1,0,column),Ip);

   Ip = -tval - Vval;
   O.insert(shape(row,3,2,column),Ip);

}

/**
 * insert complementary operator for a_up: behaves as a creator of up
 */
void insert_triple_anni_up_last(QSDArray<4> &O,int row,int column,double Vval){

   DArray<4> Ip(1,1,1,1);

   Ip = Vval;
   O.insert(shape(row,3,2,column),Ip);

}

/**
 * insert complementary operator for a_down: behaves as a creator of down 
 */
void insert_triple_anni_down_first(QSDArray<4> &O,int row,int column,double tval,double Vval){

   DArray<4> Ip(1,1,1,1);

   Ip = tval;
   O.insert(shape(row,2,0,column),Ip);

   Ip = tval + Vval;
   O.insert(shape(row,3,1,column),Ip);

}

/**
 * insert complementary operator for a_down: behaves as a creator of down 
 */
void insert_triple_anni_down_last(QSDArray<4> &O,int row,int column,double Vval){

   DArray<4> Ip(1,1,1,1);

   Ip = Vval;
   O.insert(shape(row,3,1,column),Ip);

}

/**
 * insert complementary operator for a^+_down: behaves as an annihilator of down 
 */
void insert_triple_crea_down_first(QSDArray<4> &O,int row,int column,double tval,double Vval){

   DArray<4> Ip(1,1,1,1);

   Ip = tval;
   O.insert(shape(row,0,2,column),Ip);

   Ip = tval + Vval;
   O.insert(shape(row,1,3,column),Ip);

}

/**
 * insert complementary operator for a^+_down: behaves as an annihilator of down 
 */
void insert_triple_crea_down_last(QSDArray<4> &O,int row,int column,double Vval){

   DArray<4> Ip(1,1,1,1);

   Ip = Vval;
   O.insert(shape(row,1,3,column),Ip);

}

/**
 * insert complementary operator for a^+_up: behaves as an annihilator of up
 */
void insert_triple_crea_up_first(QSDArray<4> &O,int row,int column,double tval,double Vval){

   DArray<4> Ip(1,1,1,1);

   Ip = tval;
   O.insert(shape(row,0,1,column),Ip);

   Ip = -tval - Vval;
   O.insert(shape(row,2,3,column),Ip);

}

/**
 * insert complementary operator for a^+_up: behaves as an annihilator of up
 */
void insert_triple_crea_up_last(QSDArray<4> &O,int row,int column,double Vval){

   DArray<4> Ip(1,1,1,1);

   Ip = Vval;
   O.insert(shape(row,2,3,column),Ip);

}

/**
 * insert local term: t(i,i) + V(i,i,i,i)
 */
void insert_local(QSDArray<4> &O,int row,int column,double tval,double Vval){

   DArray<4> Ip(1,1,1,1);

   Ip = tval;

   O.insert(shape(row,1,1,column),Ip);
   O.insert(shape(row,2,2,column),Ip);

   Ip = 2*tval + Vval;

   O.insert(shape(row,3,3,column),Ip);

}

/**
 * insert pair val[0] a^+_up a_up (-1)^n_down + val[1] a^+_down a_down (-1)^n_up
 */
void insert_pair_s(QSDArray<4> &O,int row,int column,const std::vector<double> &val){

   DArray<4> Ip(1,1,1,1);

   //up up
   Ip = val[0];
   O.insert(shape(row,1,1,column),Ip);

   //down down
   Ip = val[1];
   O.insert(shape(row,2,2,column),Ip);

   //both
   Ip = -val[0] - val[1];
   O.insert(shape(row,3,3,column),Ip);

}

/**
 * insert pair val[0] a^+_up a_up (-1)^n_down + val[1] a^+_down a_down (-1)^n_up
 */
void insert_pair(QSDArray<4> &O,int row,int column,const std::vector<double> &val){

   DArray<4> Ip(1,1,1,1);

   //up up
   Ip = val[0];
   O.insert(shape(row,1,1,column),Ip);

   //down down
   Ip = val[1];
   O.insert(shape(row,2,2,column),Ip);

   //both
   Ip = val[0] + val[1];
   O.insert(shape(row,3,3,column),Ip);

}

/**
 * insert general local term for one body operator
 */
void insert_local_ob(QSDArray<4> &O,int row,int column,double val){

   DArray<4> Tloc(1,1,1,1);
   Tloc = val;

   O.insert(shape(row,1,1,column),Tloc);
   O.insert(shape(row,2,2,column),Tloc);

   Tloc = 2.0*val;
   O.insert(shape(row,3,3,column),Tloc);

}

/**
 * fill the T2 array with the mp2 ansatz
 */
void fill_mp2(DArray<4> &T,const DArray<4> &V,const std::vector<double> &e){

   int no = T.shape(0);
   int nv = T.shape(2);
   int L = no + nv;

   for(int i = 0;i < no;++i)
      for(int j = 0;j < no;++j)
         for(int a = no;a < L;++a)
            for(int b = no;b < L;++b)
               T(i,j,a - no,b - no) += V(i,j,a,b) / ( e[i] + e[j] - e[a] - e[b] );

}

/**
 * given a list with all one-dimensional unmerged quantumnumbers q_in, return unique quantumnumbers in q_merged
 */
template<class Q>
void get_merged_index(const Qshapes<Q> &q_in,Qshapes<Q> &q_merged,std::vector< std::vector<int> > &ind_merged,std::vector< std::vector<int> > &inverse){

   q_merged.clear();
   ind_merged.clear();
   inverse.clear();

   inverse.resize(q_in.size());

   for(int i = 0;i < q_in.size();++i)
      inverse[i].resize(2);

   std::vector<int> v;

   q_merged.push_back(q_in[0]);

   v.push_back(0);
   ind_merged.push_back(v);

   v.clear();
   
   //row -> block and index
   inverse[0][0] = 0;//block = 0
   inverse[0][1] = 0;//index = 0

   for(int row = 1;row < q_in.size();++row){

      int block = is_in(q_in[row],q_merged);

      if(block == -1){//new quantumnumber

         q_merged.push_back(q_in[row]);

         v.push_back(row);
         ind_merged.push_back(v);

         v.clear();

         //row -> block and index
         inverse[row][0] = q_merged.size() - 1;//block = size of merged list
         inverse[row][1] = 0;//index = 0

      }
      else{//quantumnumber present in block

         ind_merged[block].push_back(row);

         inverse[row][0] = block;
         inverse[row][1] = ind_merged[block].size() - 1;

      }

   }

}

/**
 * if Quantumnumber qn is in the list, return the index of that quantumnumber in the list, if not, return -1
 */
template<class Q>
int is_in(const Q &qn,const Qshapes<Q> &qlist){

   for(int i = 0;i < qlist.size();++i)
      if(qn == qlist[i])
         return i;

   return -1;

}

/**
 * @return MPO object of length L containing a general one-body operator \sum_ab t_ab a^+_as a_bs
 */
template<class Q>
MPO<Q> qcham(const DArray<2> &t,const DArray<4> &V,bool merge){

   int L = t.shape(0);//number of occupied orbitals

   MPO<Q> mpo(L);

   Qshapes<Q> qp;
   physical(qp);

   Qshapes<Q> qz; // 0 quantum number
   qz.push_back(Q::zero());

   Qshapes<Q> qo;
   Qshapes<Q> qi;

   //first make the incoming and outgoing states:
   std::vector< Ostate > istates;

   //identity only incoming
   Ostate state;
   state.push_id();
   istates.push_back(state);
   state.clear();

   std::vector< Ostate > ostates;

   //identity
   state.push_id();
   ostates.push_back(state);
   state.clear();

   qo.push_back(Q::zero());//I

   //singles

   //a^+_up
   state.push_crea_up(0);
   ostates.push_back(state);
   state.clear();

   qo.push_back(Q(-1,0));

   //a^+_down 
   state.push_crea_down(0);
   ostates.push_back(state);
   state.clear();

   qo.push_back(Q(0,-1));

   //a_up 
   state.push_anni_up(0);
   ostates.push_back(state);
   state.clear();

   qo.push_back(Q(1,0));

   //a_down
   state.push_anni_down(0);
   ostates.push_back(state);
   state.clear();

   qo.push_back(Q(0,1));

   //doubles:

   //a^+_up a^+_down
   state.push_crea_down(0);
   state.push_crea_up(0);
   ostates.push_back(state);
   state.clear();

   qo.push_back(Q(-1,-1));

   //doubles: a^+_up a_up
   state.push_anni_up(0);
   state.push_crea_up(0);
   ostates.push_back(state);
   state.clear();

   qo.push_back(Q::zero());

   //doubles: a^+_up a_down
   state.push_anni_down(0);
   state.push_crea_up(0);
   ostates.push_back(state);
   state.clear();

   qo.push_back(Q(-1,1));

   //doubles: a^+_down a_up
   state.push_anni_up(0);
   state.push_crea_down(0);
   ostates.push_back(state);
   state.clear();

   qo.push_back(Q(1,-1));

   //doubles: a^+_down a_down
   state.push_anni_down(0);
   state.push_crea_down(0);
   ostates.push_back(state);
   state.clear();

   qo.push_back(Q::zero());

   //doubles: a_up a_down
   state.push_anni_up(0);
   state.push_anni_down(0);
   ostates.push_back(state);
   state.clear();

   qo.push_back(Q(1,1));

   //complementary operators: triples: they have the state signature of the operator they are going to, but the opposite quantumnumber
   for(int j = 1;j < L;++j){

      state.push_crea_up(j);
      ostates.push_back(state);
      state.clear();

      qo.push_back(Q(1,0));

      state.push_crea_down(j);
      ostates.push_back(state);
      state.clear();

      qo.push_back(Q(0,1));

      state.push_anni_up(j);
      ostates.push_back(state);
      state.clear();

      qo.push_back(Q(-1,0));

      state.push_anni_down(j);
      ostates.push_back(state);
      state.clear();

      qo.push_back(Q(0,-1));

   }

   //finally the local term:
   state.push_id();
   ostates.push_back(state);
   state.clear();

   qo.push_back(Q::zero());

   //now fill
   mpo[0].resize(Q::zero(),make_array(qz,qp,-qp,qo));

   int row = 0;
   int column = 0;

   //insert id
   insert_id(mpo[0],0,0);

   //insert singles
   insert_crea_up_s(mpo[0],0,1,1.0);
   insert_crea_down(mpo[0],0,2,1.0);
   insert_anni_up_s(mpo[0],0,3,1.0);
   insert_anni_down(mpo[0],0,4,1.0);

   column = 5;

   //insert doubles
   insert_crea_up_crea_down(mpo[0],0,5,1.0);
   insert_crea_up_anni_up(mpo[0],0,6,1.0);
   insert_crea_up_anni_down(mpo[0],0,7,1.0);
   insert_crea_down_anni_up(mpo[0],0,8,1.0);
   insert_crea_down_anni_down(mpo[0],0,9,1.0);
   insert_anni_down_anni_up(mpo[0],0,10,1.0);

   column = 11;

   //insert triples: complementary operator
   for(int j = 1;j < L;++j){

      insert_triple_crea_up_first(mpo[0],0,column,t(0,j),V(0,0,j,0));column++;
      insert_triple_crea_down_first(mpo[0],0,column,t(0,j),V(0,0,j,0));column++;
      insert_triple_anni_up_first(mpo[0],0,column,t(0,j),V(0,0,j,0));column++;
      insert_triple_anni_down_first(mpo[0],0,column,t(0,j),V(0,0,j,0));column++;

   }

   //last term:
   insert_local(mpo[0],0,column,t(0,0),V(0,0,0,0));

   istates = ostates;
   qi = qo;

   //middle tensors
   for(int i = 1;i < L/2;++i){

      //first construct the ingoing and outgoing states
      ostates.clear();
      qo.clear();

      //identity
      state.push_id();
      ostates.push_back(state);
      state.clear();

      qo.push_back(Q::zero());//I

      //singles

      //a^+_up
      state.push_crea_up(i);
      ostates.push_back(state);
      state.clear();

      qo.push_back(Q(-1,0));

      //a^+_down 
      state.push_crea_down(i);
      ostates.push_back(state);
      state.clear();

      qo.push_back(Q(0,-1));

      //a_up 
      state.push_anni_up(i);
      ostates.push_back(state);
      state.clear();

      qo.push_back(Q(1,0));

      //a_down
      state.push_anni_down(i);
      ostates.push_back(state);
      state.clear();

      qo.push_back(Q(0,1));

      //copy the singles from previous tensor
      row = 1;

      while(istates[row].size() == 1){

         ostates.push_back(istates[row]);
         qo.push_back(qi[row]);
         ++row;

      }

      //doubles:

      //a^+_up a^+_down
      state.push_crea_down(i);
      state.push_crea_up(i);
      ostates.push_back(state);
      state.clear();

      qo.push_back(Q(-1,-1));

      //doubles: a^+_up a_up
      state.push_anni_up(i);
      state.push_crea_up(i);
      ostates.push_back(state);
      state.clear();

      qo.push_back(Q::zero());

      //doubles: a^+_up a_down
      state.push_anni_down(i);
      state.push_crea_up(i);
      ostates.push_back(state);
      state.clear();

      qo.push_back(Q(-1,1));

      //doubles: a^+_down a_up
      state.push_anni_up(i);
      state.push_crea_down(i);
      ostates.push_back(state);
      state.clear();

      qo.push_back(Q(1,-1));

      //doubles: a^+_down a_down
      state.push_anni_down(i);
      state.push_crea_down(i);
      ostates.push_back(state);
      state.clear();

      qo.push_back(Q::zero());

      //doubles: a_down a_up 
      state.push_anni_up(i);
      state.push_anni_down(i);
      ostates.push_back(state);
      state.clear();

      qo.push_back(Q(1,1));

      //make new doubles by adding on single f^'s and f's
      row = 1;

      while(istates[row].size() == 1){//create up

         state.push_crea_up(i);
         state.insert(state.end(),istates[row].begin(),istates[row].end());
         ostates.push_back(state);
         state.clear();

         Quantum tmp = qi[row];
         tmp.crea_up();
         qo.push_back(tmp);

         ++row;

      }

      row = 1;

      while(istates[row].size() == 1){//create down

         state.push_crea_down(i);
         state.insert(state.end(),istates[row].begin(),istates[row].end());
         ostates.push_back(state);
         state.clear();

         Quantum tmp = qi[row];
         tmp.crea_down();
         qo.push_back(tmp);

         ++row;

      }

      row = 1;

      while(istates[row].size() == 1){//anni up

         state.push_anni_up(i);
         state.insert(state.end(),istates[row].begin(),istates[row].end());
         ostates.push_back(state);
         state.clear();

         Quantum tmp = qi[row];
         tmp.anni_up();
         qo.push_back(tmp);

         ++row;

      }

      row = 1;

      while(istates[row].size() == 1){//anni down 

         state.push_anni_down(i);
         state.insert(state.end(),istates[row].begin(),istates[row].end());
         ostates.push_back(state);
         state.clear();

         Quantum tmp = qi[row];
         tmp.anni_down();
         qo.push_back(tmp);

         ++row;

      }

      //copy the doubles from previous tensor
      while(istates[row].size() == 2){

         ostates.push_back(istates[row]);
         qo.push_back(qi[row]);
         ++row;

      }

      //complementary operators: triples
      for(int j = i + 1;j < L;++j){

         state.push_crea_up(j);
         ostates.push_back(state);
         state.clear();

         qo.push_back(Q(1,0));

         state.push_crea_down(j);
         ostates.push_back(state);
         state.clear();

         qo.push_back(Q(0,1));

         state.push_anni_up(j);
         ostates.push_back(state);
         state.clear();

         qo.push_back(Q(-1,0));

         state.push_anni_down(j);
         ostates.push_back(state);
         state.clear();

         qo.push_back(Q(0,-1));

      }

      //finally the local term:
      state.push_id();
      ostates.push_back(state);
      state.clear();

      qo.push_back(Q::zero());

      //now fill
      mpo[i].resize(Q::zero(),make_array(-qi,qp,-qp,qo));

      row = 0;
      column = 0;

      //insert id
      insert_id(mpo[i],0,0);

      //insert singles
      insert_crea_up_s(mpo[i],0,1,1.0);
      insert_crea_down(mpo[i],0,2,1.0);
      insert_anni_up_s(mpo[i],0,3,1.0);
      insert_anni_down(mpo[i],0,4,1.0);

      //insert signs
      row = 1;
      column = 5;

      while(istates[row].size() == 1){

         insert_sign(mpo[i],row,column);
         ++row;
         ++column;

      }

      //insert doubles
      insert_crea_up_crea_down(mpo[i],0,column,1.0);++column;
      insert_crea_up_anni_up(mpo[i],0,column,1.0);++column;
      insert_crea_up_anni_down(mpo[i],0,column,1.0);++column;
      insert_crea_down_anni_up(mpo[i],0,column,1.0);++column;
      insert_crea_down_anni_down(mpo[i],0,column,1.0);++column;
      insert_anni_down_anni_up(mpo[i],0,column,1.0);++column;

      //insert singles to form outgoing doubles
      row = 1;

      while(istates[row].size() == 1){//create up

         insert_crea_up(mpo[i],row,column,1.0);
         ++row;
         ++column;

      }

      row = 1;

      while(istates[row].size() == 1){//create down with sign

         insert_crea_down_s(mpo[i],row,column,1.0);
         ++row;
         ++column;

      }

      row = 1;

      while(istates[row].size() == 1){//annihilate up

         insert_anni_up(mpo[i],row,column,1.0);
         ++row;
         ++column;

      }

      row = 1;

      while(istates[row].size() == 1){//annihilate down with sign

         insert_anni_down_s(mpo[i],row,column,1.0);
         ++row;
         ++column;

      }

      //copy the doubles from previous tensor: identity
      while(istates[row].size() == 2){

         insert_id(mpo[i],row,column);
         ++row;
         ++column;

      }

      //HERE STARTS THE COMPLEMENTARY OPERATOR STUFF!
      while(column < ostates.size() - 1){

         int j = ostates[column].gsite(0);
         int sj = ostates[column].gspin(0);
         int aj = ostates[column].gact(0);

         //first row
         if(sj == 0 && aj == 0)
            insert_triple_crea_up_first(mpo[i],0,column,t(i,j),V(i,i,j,i));
         else if(sj == 1 && aj == 0)
            insert_triple_crea_down_first(mpo[i],0,column,t(i,j),V(i,i,j,i));
         else if(sj == 0 && aj == 1)
            insert_triple_anni_up_first(mpo[i],0,column,t(i,j),V(i,i,j,i));
         else
            insert_triple_anni_down_first(mpo[i],0,column,t(i,j),V(i,i,j,i));

         //singles coming in:
         row = 1;

         while(istates[row].size() == 1){

            std::vector<double> val(2);

            std::vector<int> v = Ostate::get_double_complement(i,istates[row],ostates[column],V,val);

            if(v.size() == 1){

               if(v[0] == 0)
                  insert_anni_down_anni_up(mpo[i],row,column,-val[0]);//extra minus sign!
               else if(v[0] == 1)
                  insert_crea_up_crea_down(mpo[i],row,column,val[0]);
               else if(v[0] == 2)
                  insert_crea_up_anni_down(mpo[i],row,column,val[0]);
               else
                  insert_crea_down_anni_up(mpo[i],row,column,val[0]);

            }
            else if(v.size() == 2)//with sign because in the middle!
               insert_pair_s(mpo[i],row,column,val);

            ++row;

         }

         //pairs coming in
         while(istates[row].size() == 2){

            double val;

            std::vector<int> v = Ostate::get_single_complement(i,istates[row],ostates[column],V,val);

            if(v.size() > 0){

               if(v[0] == 0 && v[1] == 0)//create spin up
                  insert_crea_up_s(mpo[i],row,column,val);
               else if(v[0] == 1 && v[1] == 0)//create spin down
                  insert_crea_down(mpo[i],row,column,val);
               else if(v[0] == 0 && v[1] == 1)//annihilate spin up
                  insert_anni_up_s(mpo[i],row,column,val);
               else//annihilate spin down
                  insert_anni_down(mpo[i],row,column,val);

            }

            ++row;

         }

         //signs: find row and column which are connected
         while(istates[row] != ostates[column])
            ++row;

         insert_sign(mpo[i],row,column);

         column++;

      }

      //last column! closing down everything!

      //first row
      insert_local(mpo[i],0,column,t(i,i),V(i,i,i,i));

      //close down the singles coming in with a triplet
      row = 1;

      while(istates[row].size() == 1){

         //incoming operator
         int k = istates[row].gsite(0);
         int sk = istates[row].gspin(0);
         int ak = istates[row].gact(0);

         if(sk == 0 && ak == 0)//create up coming in
            insert_triple_crea_up_last(mpo[i],row,column,V(k,i,i,i));
         else if(sk == 1 && ak == 0)//create down coming in
            insert_triple_crea_down_last(mpo[i],row,column,-V(i,k,i,i));
         else if(sk == 0 && ak == 1)//annihilate up coming in
            insert_triple_anni_up_last(mpo[i],row,column,V(i,i,k,i));
         else //annihilate down coming in
            insert_triple_anni_down_last(mpo[i],row,column,-V(i,i,i,k));

         ++row;

      }

      //close down the doubles coming in with a pair
      while(istates[row].size() == 2){

         std::vector<double> val(2);

         std::vector<int> v = Ostate::get_closing_pair_in(i,istates[row],V,val);

         if(v.size() == 1){

            if(v[0] == 0)
               insert_anni_down_anni_up(mpo[i],row,column,-val[0]);//extra minus sign!
            else if(v[0] == 1)
               insert_crea_up_crea_down(mpo[i],row,column,val[0]);
            else if(v[0] == 2)
               insert_crea_up_anni_down(mpo[i],row,column,val[0]);
            else
               insert_crea_down_anni_up(mpo[i],row,column,val[0]);

         }
         else if(v.size() == 2)
            insert_pair(mpo[i],row,column,val);

         ++row;

      }

      //close down the complementary operators of this site
      //basically the first 4 incoming states should be closed down
      insert_crea_up(mpo[i],row,column,1.0);
      ++row;

      insert_crea_down_s(mpo[i],row,column,1.0);
      ++row;

      insert_anni_up(mpo[i],row,column,1.0);
      ++row;

      insert_anni_down_s(mpo[i],row,column,1.0);
      ++row;

      //finally insert the identity on the lower right element
      insert_id(mpo[i],istates.size() - 1,ostates.size() - 1);

      istates = ostates;
      qi = qo;

   }

   //now switch around incoming and outgoing pairs to reduce the dimension of the MPO
   ostates.clear();
   qo.clear();

   //first row is closed
   state.push_id();
   ostates.push_back(state);
   state.clear();

   qo.push_back(Q::zero());//I

   //outgoing singles
   for(int j = L/2 + 1;j < L;++j){

      state.push_crea_up(j);
      ostates.push_back(state);
      state.clear();

      qo.push_back(Q(1,0));

      state.push_crea_down(j);
      ostates.push_back(state);
      state.clear();

      qo.push_back(Q(0,1));

      state.push_anni_up(j);
      ostates.push_back(state);
      state.clear();

      qo.push_back(Q(-1,0));

      state.push_anni_down(j);
      ostates.push_back(state);
      state.clear();

      qo.push_back(Q(0,-1));

   }

   //outgoing pairs: inverse quantum numbers!
   Ostate istate;

   for(int i = L/2 + 1;i < L;++i){

      //first crea up
      istate.push_crea_up(i);

      //first on-site terms

      //a^+_up a^+_down
      state = istate;
      state.push_crea_down(i);
      ostates.push_back(state);
      state.clear();

      qo.push_back(Q(1,1));

      //a^+_up a_up
      state = istate;
      state.push_anni_up(i);
      ostates.push_back(state);
      state.clear();

      qo.push_back(Q::zero());

      //a^+_up a_down
      state = istate;
      state.push_anni_down(i);
      ostates.push_back(state);
      state.clear();

      qo.push_back(Q(1,-1));

      //two-site terms
      for(int j = i + 1;j < L;++j){

         //a^+_up a^+_up
         state = istate;
         state.push_crea_up(j);
         ostates.push_back(state);
         state.clear();

         qo.push_back(Q(2,0));

         //a^+_up a^+_down
         state = istate;
         state.push_crea_down(j);
         ostates.push_back(state);
         state.clear();

         qo.push_back(Q(1,1));

         //a^+_up a_down
         state = istate;
         state.push_anni_down(j);
         ostates.push_back(state);
         state.clear();

         qo.push_back(Q(1,-1));

         //a^+_up a_up
         state = istate;
         state.push_anni_up(j);
         ostates.push_back(state);
         state.clear();

         qo.push_back(Q(0,0));

      }

      istate.clear();

      //first crea down
      istate.push_crea_down(i);

      //a^+_down a_up
      state = istate;
      state.push_anni_up(i);
      ostates.push_back(state);
      state.clear();

      qo.push_back(Q(-1,1));

      //a^+_down a_down
      state = istate;
      state.push_anni_down(i);
      ostates.push_back(state);
      state.clear();

      qo.push_back(Q::zero());

      for(int j = i + 1;j < L;++j){

         //a^+_down a^+_up
         state = istate;
         state.push_crea_up(j);
         ostates.push_back(state);
         state.clear();

         qo.push_back(Q(1,1));

         //a^+_down a^+_down
         state = istate;
         state.push_crea_down(j);
         ostates.push_back(state);
         state.clear();

         qo.push_back(Q(0,2));

         //a^+_down a_down
         state = istate;
         state.push_anni_down(j);
         ostates.push_back(state);
         state.clear();

         qo.push_back(Q(0,0));

         //a^+_down a_up
         state = istate;
         state.push_anni_up(j);
         ostates.push_back(state);
         state.clear();

         qo.push_back(Q(-1,1));

      }

      istate.clear();

      //first anni down
      istate.push_anni_down(i);

      //a_down a_up 
      state = istate;
      state.push_anni_up(i);
      ostates.push_back(state);
      state.clear();

      qo.push_back(Q(-1,-1));

      for(int j = i + 1;j < L;++j){

         //a_down a^+_up
         state = istate;
         state.push_crea_up(j);
         ostates.push_back(state);
         state.clear();

         qo.push_back(Q(1,-1));

         //a_down a^+_down
         state = istate;
         state.push_crea_down(j);
         ostates.push_back(state);
         state.clear();

         qo.push_back(Q(0,0));

         //a_down a_down
         state = istate;
         state.push_anni_down(j);
         ostates.push_back(state);
         state.clear();

         qo.push_back(Q(0,-2));

         //a_down a_up
         state = istate;
         state.push_anni_up(j);
         ostates.push_back(state);
         state.clear();

         qo.push_back(Q(-1,-1));

      }
      
      istate.clear();

      //first is anni up 
      istate.push_anni_up(i);

      for(int j = i + 1;j < L;++j){

         //a_up a^+_up
         state = istate;
         state.push_crea_up(j);
         ostates.push_back(state);
         state.clear();

         qo.push_back(Q(0,0));

         //a_up a^+_down
         state = istate;
         state.push_crea_down(j);
         ostates.push_back(state);
         state.clear();

         qo.push_back(Q(-1,1));

         //a_up a_down
         state = istate;
         state.push_anni_down(j);
         ostates.push_back(state);
         state.clear();

         qo.push_back(Q(-1,-1));

         //a_up a_up
         state = istate;
         state.push_anni_up(j);
         ostates.push_back(state);
         state.clear();

         qo.push_back(Q(-2,0));

      }

      istate.clear();

   }

   //incoming singles
   for(int i = 0;i < L/2 + 1;++i){

      //a^+_up
      state.push_crea_up(i);
      ostates.push_back(state);
      state.clear();

      qo.push_back(Q(-1,0));

      //a^+_down
      state.push_crea_down(i);
      ostates.push_back(state);
      state.clear();

      qo.push_back(Q(0,-1));

      //a_down
      state.push_anni_down(i);
      ostates.push_back(state);
      state.clear();

      qo.push_back(Q(0,1));

      //a_up
      state.push_anni_up(i);
      ostates.push_back(state);
      state.clear();

      qo.push_back(Q(1,0));

   }

   //finally last column the identity
   state.push_id();
   ostates.push_back(state);
   state.clear();

   qo.push_back(Q::zero());

   //now fill
   mpo[L/2].resize(Q::zero(),make_array(-qi,qp,-qp,qo));

   //first row completely closed
   insert_local(mpo[L/2],0,0,t(L/2,L/2),V(L/2,L/2,L/2,L/2));

   //insert complementary triples:
   row = 1;

   while(istates[row].size() == 1){

      //incoming operator
      int k = istates[row].gsite(0);
      int sk = istates[row].gspin(0);
      int ak = istates[row].gact(0);

      if(sk == 0 && ak == 0)//create up coming in
         insert_triple_crea_up_last(mpo[L/2],row,0,V(k,L/2,L/2,L/2));
      else if(sk == 1 && ak == 0)//create down coming in
         insert_triple_crea_down_last(mpo[L/2],row,0,-V(L/2,k,L/2,L/2));
      else if(sk == 0 && ak == 1)//annihilate up coming in
         insert_triple_anni_up_last(mpo[L/2],row,0,V(L/2,L/2,k,L/2));
      else //annihilate down coming in
         insert_triple_anni_down_last(mpo[L/2],row,0,-V(L/2,L/2,L/2,k));

      ++row;

   }

   //close down the doubles coming in with a pair
   while(istates[row].size() == 2){

      std::vector<double> val(2);

      std::vector<int> v = Ostate::get_closing_pair_in(L/2,istates[row],V,val);

      if(v.size() == 1){

         if(v[0] == 0)
            insert_anni_down_anni_up(mpo[L/2],row,0,-val[0]);//extra minus sign!
         else if(v[0] == 1)
            insert_crea_up_crea_down(mpo[L/2],row,0,val[0]);
         else if(v[0] == 2)
            insert_crea_up_anni_down(mpo[L/2],row,0,val[0]);
         else
            insert_crea_down_anni_up(mpo[L/2],row,0,val[0]);

      }
      else if(v.size() == 2)
         insert_pair(mpo[L/2],row,0,val);

      ++row;

   }

   //close down the complementary triple operators of this site
   //basically the first 4 incoming states should be closed down
   insert_crea_up(mpo[L/2],row,0,1.0);
   ++row;

   insert_crea_down_s(mpo[L/2],row,0,1.0);
   ++row;

   insert_anni_up(mpo[L/2],row,0,1.0);
   ++row;

   insert_anni_down_s(mpo[L/2],row,0,1.0);
   ++row;

   //finally insert the identity for the already closed terms coming in
   insert_id(mpo[L/2],istates.size() - 1,0);

   //next columns: outgoing singles
   int col = 1;

   while(ostates[col].size() == 1){

      int j = ostates[col].gsite(0);
      int sj = ostates[col].gspin(0);
      int aj = ostates[col].gact(0);

      //first row
      if(sj == 0 && aj == 0)
         insert_triple_crea_up_first(mpo[L/2],0,col,t(L/2,j),V(L/2,L/2,j,L/2));
      else if(sj == 1 && aj == 0)
         insert_triple_crea_down_first(mpo[L/2],0,col,t(L/2,j),V(L/2,L/2,j,L/2));
      else if(sj == 0 && aj == 1)
         insert_triple_anni_up_first(mpo[L/2],0,col,t(L/2,j),V(L/2,L/2,j,L/2));
      else
         insert_triple_anni_down_first(mpo[L/2],0,col,t(L/2,j),V(L/2,L/2,j,L/2));

      //singles coming in:
      row = 1;

      while(istates[row].size() == 1){

         std::vector<double> val(2);

         std::vector<int> v = Ostate::get_double_complement(L/2,istates[row],ostates[col],V,val);

         if(v.size() == 1){

            if(v[0] == 0)
               insert_anni_down_anni_up(mpo[L/2],row,col,-val[0]);//extra minus sign!
            else if(v[0] == 1)
               insert_crea_up_crea_down(mpo[L/2],row,col,val[0]);
            else if(v[0] == 2)
               insert_crea_up_anni_down(mpo[L/2],row,col,val[0]);
            else
               insert_crea_down_anni_up(mpo[L/2],row,col,val[0]);

         }
         else if(v.size() == 2)//with sign because in the middle!
            insert_pair_s(mpo[L/2],row,col,val);

         ++row;

      }

      //pairs coming in
      while(istates[row].size() == 2){

         double val;

         std::vector<int> v = Ostate::get_single_complement(L/2,istates[row],ostates[col],V,val);

         if(v.size() > 0){

            if(v[0] == 0 && v[1] == 0)//create spin up
               insert_crea_up_s(mpo[L/2],row,col,val);
            else if(v[0] == 1 && v[1] == 0)//create spin down
               insert_crea_down(mpo[L/2],row,col,val);
            else if(v[0] == 0 && v[1] == 1)//annihilate spin up
               insert_anni_up_s(mpo[L/2],row,col,val);
            else//annihilate spin down
               insert_anni_down(mpo[L/2],row,col,val);

         }

         ++row;

      }

      //signs: find row and column which are connected
      while(istates[row] != ostates[col])
         ++row;

      insert_sign(mpo[L/2],row,col);

      col++;

   }

   //now outgoing pairs!
   while(ostates[col].size() == 2){

      int k = ostates[col].gsite(0);
      int sk = ostates[col].gspin(0);
      int ak = ostates[col].gact(0);

      int l = ostates[col].gsite(1);
      int sl = ostates[col].gspin(1);
      int al = ostates[col].gact(1);

      //first row: id coming in, insert pairs
      std::vector<double> val(2);

      std::vector<int> v = Ostate::get_closing_pair_out(L/2,ostates[col],V,val);

      if(v.size() == 1){

         if(v[0] == 0)
            insert_anni_down_anni_up(mpo[L/2],0,col,-val[0]);//extra minus sign!
         else if(v[0] == 1)
            insert_crea_up_crea_down(mpo[L/2],0,col,val[0]);
         else if(v[0] == 2)
            insert_crea_up_anni_down(mpo[L/2],0,col,val[0]);
         else
            insert_crea_down_anni_up(mpo[L/2],0,col,val[0]);

      }
      else if(v.size() == 2)
         insert_pair(mpo[L/2],0,col,val);

      //singles coming in: insert single to match outgoing pair
      row = 1;

      while(istates[row].size() == 1){

         double val;

         std::vector<int> v = Ostate::get_single_complement(L/2,ostates[col],istates[row],V,val);

         if(v.size() > 0){

            if(v[0] == 0 && v[1] == 0)//create spin up
               insert_crea_up(mpo[L/2],row,col,val);
            else if(v[0] == 1 && v[1] == 0)//create spin down
               insert_crea_down_s(mpo[L/2],row,col,val);
            else if(v[0] == 0 && v[1] == 1)//annihilate spin up
               insert_anni_up(mpo[L/2],row,col,val);
            else//annihilate spin down
               insert_anni_down_s(mpo[L/2],row,col,val);

         }

         ++row;

      }

      //pairs coming in: convert incoming to outgoing pairs

      while(istates[row].size() == 2){

         double val;

         int op = Ostate::transfer_pair_in_pair_out(istates[row],ostates[col],V,val);

         if(op == 1)
            insert_id(mpo[L/2],row,col,val);

         ++row;

      }

      ++col;

   }

   //then the outgoing "incoming" singles
   while(col < ostates.size() - 1){

      int ci = ostates[col].gsite(0);
      int cs = ostates[col].gspin(0);
      int ca = ostates[col].gact(0);

      //insert single operators on the first row
      if(ci == L/2){

         if(cs == 0 && ca == 0)
            insert_crea_up_s(mpo[L/2],0,col,1.0);
         else if(cs == 1 && ca == 0)
            insert_crea_down(mpo[L/2],0,col,1.0);
         else if(cs == 0 && ca == 1)
            insert_anni_up_s(mpo[L/2],0,col,1.0);
         else
            insert_anni_down(mpo[L/2],0,col,1.0);

      }

      row = 1;

      while(istates[row].size() == 1){

         int ri = istates[row].gsite(0);
         int rs = istates[row].gspin(0);
         int ra = istates[row].gact(0);

         if(ri == ci && rs == cs && ra == ca)
            insert_sign(mpo[L/2],row,col);

         ++row;

      }

      ++col;

   }

   //finally insert id for incoming id
   insert_id(mpo[L/2],0,ostates.size() - 1);

   istates = ostates;
   qi = qo;

   //the rest of the blocks: L/2 + 1 --> L - 1
   for(int i = L/2 + 1;i < L - 1;++i){

      ostates.clear();
      qo.clear();

      //first row is closed
      state.push_id();
      ostates.push_back(state);
      state.clear();

      qo.push_back(Q::zero());//I

      //outgoing singles
      for(int j = i + 1;j < L;++j){

         state.push_crea_up(j);
         ostates.push_back(state);
         state.clear();

         qo.push_back(Q(1,0));

         state.push_crea_down(j);
         ostates.push_back(state);
         state.clear();

         qo.push_back(Q(0,1));

         state.push_anni_up(j);
         ostates.push_back(state);
         state.clear();

         qo.push_back(Q(-1,0));

         state.push_anni_down(j);
         ostates.push_back(state);
         state.clear();

         qo.push_back(Q(0,-1));

      }

      //outgoing pairs: inverse quantum numbers!
      Ostate istate;

      for(int j = i + 1;j < L;++j){

         //first crea up
         istate.push_crea_up(j);

         //first on-site terms

         //a^+_up a^+_down
         state = istate;
         state.push_crea_down(j);
         ostates.push_back(state);
         state.clear();

         qo.push_back(Q(1,1));

         //a^+_up a_up
         state = istate;
         state.push_anni_up(j);
         ostates.push_back(state);
         state.clear();

         qo.push_back(Q::zero());

         //a^+_up a_down
         state = istate;
         state.push_anni_down(j);
         ostates.push_back(state);
         state.clear();

         qo.push_back(Q(1,-1));

         //two-site terms
         for(int k = j + 1;k < L;++k){

            //a^+_up a^+_up
            state = istate;
            state.push_crea_up(k);
            ostates.push_back(state);
            state.clear();

            qo.push_back(Q(2,0));

            //a^+_up a^+_down
            state = istate;
            state.push_crea_down(k);
            ostates.push_back(state);
            state.clear();

            qo.push_back(Q(1,1));

            //a^+_up a_down
            state = istate;
            state.push_anni_down(k);
            ostates.push_back(state);
            state.clear();

            qo.push_back(Q(1,-1));

            //a^+_up a_up
            state = istate;
            state.push_anni_up(k);
            ostates.push_back(state);
            state.clear();

            qo.push_back(Q(0,0));

         }

         istate.clear();

         //first crea down
         istate.push_crea_down(j);

         //a^+_down a_up
         state = istate;
         state.push_anni_up(j);
         ostates.push_back(state);
         state.clear();

         qo.push_back(Q(-1,1));

         //a^+_down a_down
         state = istate;
         state.push_anni_down(j);
         ostates.push_back(state);
         state.clear();

         qo.push_back(Q::zero());

         for(int k = j + 1;k < L;++k){

            //a^+_down a^+_up
            state = istate;
            state.push_crea_up(k);
            ostates.push_back(state);
            state.clear();

            qo.push_back(Q(1,1));

            //a^+_down a^+_down
            state = istate;
            state.push_crea_down(k);
            ostates.push_back(state);
            state.clear();

            qo.push_back(Q(0,2));

            //a^+_down a_down
            state = istate;
            state.push_anni_down(k);
            ostates.push_back(state);
            state.clear();

            qo.push_back(Q(0,0));

            //a^+_down a_up
            state = istate;
            state.push_anni_up(k);
            ostates.push_back(state);
            state.clear();

            qo.push_back(Q(-1,1));

         }

         istate.clear();

         //first anni down
         istate.push_anni_down(j);

         //a_down a_up 
         state = istate;
         state.push_anni_up(j);
         ostates.push_back(state);
         state.clear();

         qo.push_back(Q(-1,-1));

         for(int k = j + 1;k < L;++k){

            //a_down a^+_up
            state = istate;
            state.push_crea_up(k);
            ostates.push_back(state);
            state.clear();

            qo.push_back(Q(1,-1));

            //a_down a^+_down
            state = istate;
            state.push_crea_down(k);
            ostates.push_back(state);
            state.clear();

            qo.push_back(Q(0,0));

            //a_down a_down
            state = istate;
            state.push_anni_down(k);
            ostates.push_back(state);
            state.clear();

            qo.push_back(Q(0,-2));

            //a_down a_up
            state = istate;
            state.push_anni_up(k);
            ostates.push_back(state);
            state.clear();

            qo.push_back(Q(-1,-1));

         }

         istate.clear();

         //first is anni up 
         istate.push_anni_up(j);

         for(int k = j + 1;k < L;++k){

            //a_up a^+_up
            state = istate;
            state.push_crea_up(k);
            ostates.push_back(state);
            state.clear();

            qo.push_back(Q(0,0));

            //a_up a^+_down
            state = istate;
            state.push_crea_down(k);
            ostates.push_back(state);
            state.clear();

            qo.push_back(Q(-1,1));

            //a_up a_down
            state = istate;
            state.push_anni_down(k);
            ostates.push_back(state);
            state.clear();

            qo.push_back(Q(-1,-1));

            //a_up a_up
            state = istate;
            state.push_anni_up(k);
            ostates.push_back(state);
            state.clear();

            qo.push_back(Q(-2,0));

         }

         istate.clear();

      }

      //incoming singles
      for(int j = 0;j <= i;++j){

         //a^+_up
         state.push_crea_up(j);
         ostates.push_back(state);
         state.clear();

         qo.push_back(Q(-1,0));

         //a^+_down
         state.push_crea_down(j);
         ostates.push_back(state);
         state.clear();

         qo.push_back(Q(0,-1));

         //a_down
         state.push_anni_down(j);
         ostates.push_back(state);
         state.clear();

         qo.push_back(Q(0,1));

         //a_up
         state.push_anni_up(j);
         ostates.push_back(state);
         state.clear();

         qo.push_back(Q(1,0));

      }

      //finally last column the identity
      state.push_id();
      ostates.push_back(state);
      state.clear();

      qo.push_back(Q::zero());

      //now fill
      mpo[i].resize(Q::zero(),make_array(-qi,qp,-qp,qo));

      col = 0;

      //first column is closed
      insert_id(mpo[i],0,0);

      //for incoming complementary triples, close down with single
      row = 1;

      while(istates[row].size() == 1){

         int ri = istates[row].gsite(0);
         int rs = istates[row].gspin(0);
         int ra = istates[row].gact(0);

         if(ri == i){

            if(rs == 0 && ra == 0)
               insert_crea_up(mpo[i],row,col,1.0);
            else if(rs == 1 && ra == 0)
               insert_crea_down_s(mpo[i],row,col,1.0);
            else if(rs == 0 && ra == 1)
               insert_anni_up(mpo[i],row,col,1.0);
            else
               insert_anni_down_s(mpo[i],row,col,1.0);

         }

         ++row;

      }

      //close down incoming complementary doubles with pairs
      while(istates[row].size() == 2){

         int ri_1 = istates[row].gsite(0);
         int rs_1 = istates[row].gspin(0);
         int ra_1 = istates[row].gact(0);

         int ri_2 = istates[row].gsite(1);
         int rs_2 = istates[row].gspin(1);
         int ra_2 = istates[row].gact(1);

         if(ri_1 == i && ri_2 == i){

            if(ra_1 == 0){//first crea

               if(rs_1 == 0){//first create up

                  if(rs_2 == 1 && ra_2 == 0)
                     insert_crea_up_crea_down(mpo[i],row,col,1.0);
                  else if(rs_2 == 0 && ra_2 == 1)
                     insert_crea_up_anni_up(mpo[i],row,col,1.0);
                  else if(rs_2 == 1 && ra_2 == 1)
                     insert_crea_up_anni_down(mpo[i],row,col,1.0);

               }
               else{//first create down

                  if(rs_2 == 0 && ra_2 == 1)
                     insert_crea_down_anni_up(mpo[i],row,col,1.0);
                  else if(rs_2 == 1 && ra_2 == 1)
                     insert_crea_down_anni_down(mpo[i],row,col,1.0);

               }

            }
            else{//first anni

               if(rs_1 == 1){//first anni down

                  if(rs_2 == 0 && ra_2 == 1)//second anni up: extra minus sign!
                     insert_anni_down_anni_up(mpo[i],row,col,1.0);

               }

            }

         }

         ++row;

      }

      //incoming singles, close with triples!
      while(row < istates.size() - 1){

         //incoming operator
         int k = istates[row].gsite(0);
         int sk = istates[row].gspin(0);
         int ak = istates[row].gact(0);

         if(sk == 0 && ak == 0)//create up coming in
            insert_triple_crea_up_last(mpo[i],row,0,V(k,i,i,i));
         else if(sk == 1 && ak == 0)//create down coming in
            insert_triple_crea_down_last(mpo[i],row,0,-V(i,k,i,i));
         else if(sk == 0 && ak == 1)//annihilate up coming in
            insert_triple_anni_up_last(mpo[i],row,0,V(i,i,k,i));
         else //annihilate down coming in
            insert_triple_anni_down_last(mpo[i],row,0,-V(i,i,i,k));

         ++row;

      }

      //insert local terms
      insert_local(mpo[i],istates.size() - 1,0,t(i,i),V(i,i,i,i));

      //outgoing complementary triples
      col = 1;

      while(ostates[col].size() == 1){

         int ci = ostates[col].gsite(0);
         int cs = ostates[col].gspin(0);
         int ca = ostates[col].gact(0);

         //incoming complementary triples
         row = 1;

         while(istates[row].size() == 1){

            int ri = istates[row].gsite(0);
            int rs = istates[row].gspin(0);
            int ra = istates[row].gact(0);

            if(ci == ri && cs == rs && ca == ra)
               insert_sign(mpo[i],row,col,1.0);

            ++row;

         }

         //incoming complementary pairs: add single
         while(istates[row].size() == 2){

            int ri1 = istates[row].gsite(0);
            int rs1 = istates[row].gspin(0);
            int ra1 = istates[row].gact(0);

            int ri2 = istates[row].gsite(1);
            int rs2 = istates[row].gspin(1);
            int ra2 = istates[row].gact(1);

            if(ri1 == i){//if first site of the complementary is this site i:

               if(ri2 == ci && rs2 == cs && ra2 == ca){//and the second matches up with the outgoing single

                  if(rs1 == 0 && ra1 == 0)
                     insert_crea_up_s(mpo[i],row,col,1.0);
                  else if(rs1 == 1 && ra1 == 0)
                     insert_crea_down(mpo[i],row,col,1.0);
                  else if(rs1 == 0 && ra1 == 1)
                     insert_anni_up_s(mpo[i],row,col,1.0);
                  else
                     insert_anni_down(mpo[i],row,col,1.0);

               }

            }

            ++row;

         }

         //singles coming in:
         while(row < istates.size() - 1){

            std::vector<double> val(2);

            std::vector<int> v = Ostate::get_double_complement(i,istates[row],ostates[col],V,val);

            if(v.size() == 1){

               if(v[0] == 0)
                  insert_anni_down_anni_up(mpo[i],row,col,-val[0]);//extra minus sign!
               else if(v[0] == 1)
                  insert_crea_up_crea_down(mpo[i],row,col,val[0]);
               else if(v[0] == 2)
                  insert_crea_up_anni_down(mpo[i],row,col,val[0]);
               else
                  insert_crea_down_anni_up(mpo[i],row,col,val[0]);

            }
            else if(v.size() == 2)//with sign because in the middle!
               insert_pair_s(mpo[i],row,col,val);

            ++row;

         }

         if(cs == 0 && ca == 0)
            insert_triple_crea_up_first(mpo[i],row,col,t(i,ci),V(i,i,ci,i));
         else if(cs == 1 && ca == 0)
            insert_triple_crea_down_first(mpo[i],row,col,t(i,ci),V(i,i,ci,i));
         else if(cs == 0 && ca == 1)
            insert_triple_anni_up_first(mpo[i],row,col,t(i,ci),V(i,i,ci,i));
         else
            insert_triple_anni_down_first(mpo[i],row,col,t(i,ci),V(i,i,ci,i));

         ++col;

      }

      while(ostates[col].size() == 2){

         //now start from the bottom, with the last row
         row = istates.size() - 1;

         //last row: id coming in, insert pairs
         std::vector<double> val(2);

         std::vector<int> v = Ostate::get_closing_pair_out(i,ostates[col],V,val);

         if(v.size() == 1){

            if(v[0] == 0)
               insert_anni_down_anni_up(mpo[i],row,col,-val[0]);//extra minus sign!
            else if(v[0] == 1)
               insert_crea_up_crea_down(mpo[i],row,col,val[0]);
            else if(v[0] == 2)
               insert_crea_up_anni_down(mpo[i],row,col,val[0]);
            else
               insert_crea_down_anni_up(mpo[i],row,col,val[0]);

         }
         else if(v.size() == 2)
            insert_pair(mpo[i],row,col,val);

         row--;

         //now decrease the row: incoming singles
         while(istates[row].size() == 1){

            double val;

            std::vector<int> v = Ostate::get_single_complement(i,ostates[col],istates[row],V,val);

            if(v.size() > 0){

               if(v[0] == 0 && v[1] == 0)//create spin up
                  insert_crea_up(mpo[i],row,col,val);
               else if(v[0] == 1 && v[1] == 0)//create spin down
                  insert_crea_down_s(mpo[i],row,col,val);
               else if(v[0] == 0 && v[1] == 1)//annihilate spin up
                  insert_anni_up(mpo[i],row,col,val);
               else//annihilate spin down
                  insert_anni_down_s(mpo[i],row,col,val);

            }

            --row;

         }

         //incoming complementary pairs: put in id on the right place
         while(istates[row].size() == 2){

            if(istates[row] == ostates[col])
               insert_id(mpo[i],row,col);

            --row;

         }

         ++col;

      }

      //outgoing singles
      while(col < ostates.size() - 1){

         int ci = ostates[col].gsite(0);
         int cs = ostates[col].gspin(0);
         int ca = ostates[col].gact(0);

         row = istates.size() - 1;

         //insert correct operator:
         if(ci == i){

            if(cs == 0 && ca == 0)
               insert_crea_up_s(mpo[i],row,col,1.0);
            else if(cs == 1 && ca == 0)
               insert_crea_down(mpo[i],row,col,1.0);
            else if(cs == 0 && ca == 1)
               insert_anni_up_s(mpo[i],row,col,1.0);
            else
               insert_anni_down(mpo[i],row,col,1.0);

         }

         row--;

         //insert sign for incoming singles to transfer to outgoing singles
         while(istates[row].size() == 1){

            if(istates[row] == ostates[col])
               insert_sign(mpo[i],row,col);

            --row;

         }

         ++col;

      }

      //finally insert id in bottom right
      insert_id(mpo[i],istates.size() - 1,ostates.size() - 1);

      //transfer for following site
      qi = qo;
      istates = ostates;

   }

   //finally the last mpo:
   mpo[L - 1].resize(Q::zero(),make_array(-qi,qp,-qp,qz));

   insert_id(mpo[L - 1],0,0);

   //for incoming complementary triples, close down with single
   row = 1;

   while(istates[row].size() == 1){

      int ri = istates[row].gsite(0);
      int rs = istates[row].gspin(0);
      int ra = istates[row].gact(0);

      if(ri == L - 1){

         if(rs == 0 && ra == 0)
            insert_crea_up(mpo[L - 1],row,0,1.0);
         else if(rs == 1 && ra == 0)
            insert_crea_down_s(mpo[L - 1],row,0,1.0);
         else if(rs == 0 && ra == 1)
            insert_anni_up(mpo[L - 1],row,0,1.0);
         else
            insert_anni_down_s(mpo[L - 1],row,0,1.0);

      }

      ++row;

   }

   //close down incoming complementary doubles with pairs
   while(istates[row].size() == 2){

      int ri_1 = istates[row].gsite(0);
      int rs_1 = istates[row].gspin(0);
      int ra_1 = istates[row].gact(0);

      int ri_2 = istates[row].gsite(1);
      int rs_2 = istates[row].gspin(1);
      int ra_2 = istates[row].gact(1);

      if(ri_1 == L - 1 && ri_2 == L - 1){

         if(ra_1 == 0){//first crea

            if(rs_1 == 0){//first create up

               if(rs_2 == 1 && ra_2 == 0)
                  insert_crea_up_crea_down(mpo[L - 1],row,0,1.0);
               else if(rs_2 == 0 && ra_2 == 1)
                  insert_crea_up_anni_up(mpo[L - 1],row,0,1.0);
               else if(rs_2 == 1 && ra_2 == 1)
                  insert_crea_up_anni_down(mpo[L - 1],row,0,1.0);

            }
            else{//first create down

               if(rs_2 == 0 && ra_2 == 1)
                  insert_crea_down_anni_up(mpo[L - 1],row,0,1.0);
               else if(rs_2 == 1 && ra_2 == 1)
                  insert_crea_down_anni_down(mpo[L - 1],row,0,1.0);

            }

         }
         else{//first anni

            if(rs_1 == 1){//first anni down

               if(rs_2 == 0 && ra_2 == 1)//second anni up: extra minus sign!
                  insert_anni_down_anni_up(mpo[L - 1],row,0,1.0);

            }

         }

      }

      ++row;

   }

   //incoming singles, close with triples!
   while(row < istates.size() - 1){

      //incoming operator
      int k = istates[row].gsite(0);
      int sk = istates[row].gspin(0);
      int ak = istates[row].gact(0);

      if(sk == 0 && ak == 0)//create up coming in
         insert_triple_crea_up_last(mpo[L - 1],row,0,V(k,L - 1,L - 1,L - 1));
      else if(sk == 1 && ak == 0)//create down coming in
         insert_triple_crea_down_last(mpo[L - 1],row,0,-V(L - 1,k,L - 1,L - 1));
      else if(sk == 0 && ak == 1)//annihilate up coming in
         insert_triple_anni_up_last(mpo[L - 1],row,0,V(L - 1,L - 1,k,L - 1));
      else //annihilate down coming in
         insert_triple_anni_down_last(mpo[L - 1],row,0,-V(L - 1,L - 1,L - 1,k));

      ++row;

   }

   //insert local terms
   insert_local(mpo[L - 1],istates.size() - 1,0,t(L - 1,L - 1),V(L - 1,L - 1,L - 1,L - 1));

   if(merge == true){

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

   }

   return mpo;

}

template void physical<Quantum>(Qshapes<Quantum> &);
template MPO<Quantum> E(int,int,int,double);
template MPO<Quantum> E(int,int,int,int,int,double);
template MPO<Quantum> tpint(int,int,int,int,int,double);
template MPO<Quantum> T1(const DArray<2> &,bool);
template MPO<Quantum> T2(const DArray<4> &,bool);
template MPO<Quantum> T2_test(const DArray<4> &);
template MPO<Quantum> qcham(const DArray<2> &,const DArray<4> &,bool);
template MPO<Quantum> one_body(const DArray<2> &,bool);
template MPO<Quantum> crea_up(int L,int i);
template MPO<Quantum> crea_down(int L,int i);
template MPO<Quantum> anni_up(int L,int i);
template MPO<Quantum> anni_down(int L,int i);
template MPO<Quantum> crea_up_crea_up(int L,int i,int j,double val);
template MPO<Quantum> crea_up_crea_down(int L,int i,int j,double val);
template MPO<Quantum> crea_down_crea_up(int L,int i,int j,double val);
template MPO<Quantum> crea_down_crea_down(int L,int i,int j,double val);
template MPO<Quantum> anni_up_anni_up(int L,int i,int j,double val);
template MPO<Quantum> anni_down_anni_up(int L,int i,int j,double val);
template MPO<Quantum> anni_up_anni_down(int L,int i,int j,double val);
template MPO<Quantum> anni_down_anni_down(int L,int i,int j,double val);
template void get_merged_index(const Qshapes<Quantum> &,Qshapes<Quantum> &,std::vector< std::vector<int> > &,std::vector< std::vector<int> > &);
template int is_in(const Quantum &qn,const Qshapes<Quantum> &qlist);
