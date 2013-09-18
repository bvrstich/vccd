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
T1_2_mpo::T1_2_mpo(int no,int nv){

   this->no = no;
   this->nv = nv;

   s2ia = new int * [no*nv];

   for(int i = 0;i < no*nv;++i)
      s2ia[i] = new int [2];

   ia2s = new int * [no];

   for(int i = 0;i < no;++i)
      ia2s[i] = new int [nv];

   int s = 0;

   for(int i = 0;i < no;++i)
      for(int a = 0;a < nv;++a){

         ia2s[i][a] = s;

         s2ia[s][0] = i;
         s2ia[s][1] = a;

         ++s;

      }

   list = new std::vector< vector<int> > [no*nv];

   //now construct the list
   vector<int> mpodat(5);

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

      //first column is closed
      for(int row = 0;row < istates.size();++row){

         int ii = istates[row].gsite(0);
         int si = istates[row].gspin(0);

         if(si == 0){

            int s = ia2s[ii][0];

            mpodat[0] = no;//site
            mpodat[1] = row;//left virtual
            mpodat[2] = ;//physical outgoing
            mpodat[3] = ;//physical incoming
            mpodat[4] = 0;//right virtual

            insert_crea_up(s,row,0,t(ii,0));

         }
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

}

/**
 * copy constructor
 */
T1_2_mpo::T1_2_mpo(const T1_2_mpo &copy){

}

/**
 * destructor
 */
T1_2_mpo::~T1_2_mpo(){ }

/**
 * 
 */
