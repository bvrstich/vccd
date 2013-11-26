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
T2_2_mpo::T2_2_mpo(int no,int nv){

   this->no = no;
   this->nv = nv;

   int L = no + nv;

   //two occupied
   ij2o = new int * [no];

   for(int i = 0;i < no;++i)
      ij2o[i] = new int [no];

   int o = 0;
   std::vector<int> v(2);

   for(int i = 0;i < no;++i)
      for(int j = i;j < no;++j){

         ij2o[i][j] = o;
         ij2o[j][i] = o;

         ++o;

         v[0] = i;
         v[1] = j;

         o2ij.push_back(v);

      }

   //two virtuals
   ab2v = new int * [nv];

   for(int a = 0;a < nv;++a)
      ab2v[a] = new int [nv];

   int vi = 0;

   for(int a = 0;a < nv;++a)
      for(int b = 0;b < nv;++b){

         ab2v[a][b] = vi;

         v[0] = a;
         v[1] = b;

         v2ab.push_back(v);

         ++vi;

      }

   ov2s = new int * [no*(no + 1)/2];

   for(int o = 0;o < no*(no + 1)/2;++o)
      ov2s[o] = new int [nv*nv];

   int s = 0;

   for(int o = 0;o < no*(no + 1)/2;++o)
      for(int vi = 0;vi < nv*nv;++vi){

         ov2s[o][vi] = s;
         ++s;

         v[0] = o;
         v[1] = vi;

         s2ov.push_back(v);

      }

   list = new std::vector< std::vector<int> > [s2ov.size()];

   Qshapes<Quantum> qz;
   qz.push_back(Quantum::zero());//I

   //quantum numbers:
   Qshapes<Quantum> qo;
   qo.push_back(Quantum::zero());//I
   qo.push_back(Quantum(0,1));//a_down
   qo.push_back(Quantum(1,0));//a_up
   qo.push_back(Quantum(1,1));//a_down a_up

   Qshapes<Quantum> q_merged;
   std::vector< std::vector<int> > ind_merged;
   std::vector< std::vector<int> > inverse;

   //get merged row
   get_merged_index(qz,q_merged,ind_merged,inverse);
   merged_row.push_back(inverse);

   //get merged column
   get_merged_index(qo,q_merged,ind_merged,inverse);
   merged_col.push_back(inverse);

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

   std::vector< Ostate > istates;
   istates = ostates;

   Qshapes<Quantum> qi = qo;

   if(nv < no){

      for(int i = 1;i < nv;++i){

         ostates.clear();
         qo.clear();

         qo.push_back(Quantum::zero());//I

         ostates.push_back(istates[0]);

         qo.push_back(Quantum(0,1));//a_down
         qo.push_back(Quantum(1,0));//a_up

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

         qo.push_back(Quantum(1,1));//a_up a_down

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

         //get merged row
         get_merged_index(-qi,q_merged,ind_merged,inverse);
         merged_row.push_back(inverse);

         //get merged column
         get_merged_index(qo,q_merged,ind_merged,inverse);
         merged_col.push_back(inverse);

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

         qo.push_back(Quantum(1,1));

         for(int j = i + 1;j < L;++j){

            //up up
            state = istate;
            state.push_crea_up(j);
            ostates.push_back(state);
            state.clear();

            qo.push_back(Quantum(2,0));

            //up down
            state = istate;
            state.push_crea_down(j);
            ostates.push_back(state);
            state.clear();

            qo.push_back(Quantum(1,1));

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

            qo.push_back(Quantum(1,1));

            //down down
            state = istate;
            state.push_crea_down(j);
            ostates.push_back(state);
            state.clear();

            qo.push_back(Quantum(0,2));

         }

         istate.clear();

      }

      //then the incoming singles
      qo.push_back(Quantum(0,1));//a_down
      qo.push_back(Quantum(1,0));//a_up

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
      qo.push_back(Quantum::zero());//I
      ostates.push_back(istates[0]);

      //get merged row
      get_merged_index(-qi,q_merged,ind_merged,inverse);
      merged_row.push_back(inverse);

      //get merged column
      get_merged_index(qo,q_merged,ind_merged,inverse);
      merged_col.push_back(inverse);

      //first row
      int col = 0;

      while(ostates[col].size() == 2){

         int ai = ostates[col].gsite(0);
         int as = ostates[col].gspin(0);

         int bi = ostates[col].gsite(1);
         int bs = ostates[col].gspin(1);

         int o = ij2o[nv][nv];
         int v = ab2v[ai-no][bi-no];

         int s = ov2s[o][v];

         if(as == 0 && bs == 1)
            push_anni_down_anni_up(s,nv,0,col,1);
         else if(as == 1 && bs == 0)
            push_anni_down_anni_up(s,nv,0,col,-1);

         ++col;

      }

      //singles coming in 
      row = 1;

      while(istates[row].size() == 1){

         //to outgoing doubles
         col = 0;

         while(ostates[col].size() == 2){

            push_single_in_complement(nv,istates[row],row,ostates[col],col);

            ++col;

         }

         ++row;

      }

      //doubles coming in to doubles going out
      while(row < istates.size()){

         col = 0;

         while(ostates[col].size() == 2){

            push_double_complement(nv,istates[row],row,ostates[col],col);

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
         qo.push_back(Quantum(0,1));//a_down
         qo.push_back(Quantum(1,0));//a_up

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
         qo.push_back(Quantum::zero());//I
         ostates.push_back(istates[istates.size() - 1]);

         //get merged row
         get_merged_index(-qi,q_merged,ind_merged,inverse);
         merged_row.push_back(inverse);

         //get merged column
         get_merged_index(qo,q_merged,ind_merged,inverse);
         merged_col.push_back(inverse);

         row = 0;

         //incoming 'outgoing' doubles
         while(istates[row].size() == 2)
            ++row;

         //incoming singles
         while(row < istates.size() - 1){

            //to outgoing doubles
            col = 0;

            while(ostates[col].size() == 2){

               push_single_in_complement(i,istates[row],row,ostates[col],col);

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

            int o = ij2o[i][i];
            int v = ab2v[ai-no][bi-no];

            int s = ov2s[o][v];

            if(as == 0 && bs == 1)
               push_anni_down_anni_up(s,i,istates.size() - 1,col,1);
            else if(as == 1 && bs == 0)
               push_anni_down_anni_up(s,i,istates.size() - 1,col,-1);

            ++col;

         }

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

      //get merged row
      get_merged_index(-qi,q_merged,ind_merged,inverse);
      merged_row.push_back(inverse);

      //get merged column
      get_merged_index(qo,q_merged,ind_merged,inverse);
      merged_col.push_back(inverse);

      //incoming singles
      while(row < istates.size() - 1){

         //to outgoing doubles
         col = 0;

         while(ostates[col].size() == 2){

            push_single_in_complement(no - 1,istates[row],row,ostates[col],col);

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

         int o = ij2o[no - 1][no - 1];
         int v = ab2v[ai-no][bi-no];

         int s = ov2s[o][v];

         if(as == 0 && bs == 1)
            push_anni_down_anni_up(s,no - 1,istates.size() - 1,col,1);
         else if(as == 1 && bs == 0)
            push_anni_down_anni_up(s,no - 1,istates.size() - 1,col,-1);

         ++col;

         ++col;

      }

   }
   else{//nv >= no

      for(int i = 1;i < no - 1;++i){

         ostates.clear();
         qo.clear();

         qo.push_back(Quantum::zero());//I

         ostates.push_back(istates[0]);

         qo.push_back(Quantum(0,1));//a_down
         qo.push_back(Quantum(1,0));//a_up

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

         qo.push_back(Quantum(1,1));//a_up a_down

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

         //get merged row
         get_merged_index(-qi,q_merged,ind_merged,inverse);
         merged_row.push_back(inverse);

         //get merged column
         get_merged_index(qo,q_merged,ind_merged,inverse);
         merged_col.push_back(inverse);

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

            qo.push_back(Quantum(1,1));

            for(int j = i + 1;j < L;++j){

               //up up
               state = istate;
               state.push_crea_up(j);
               ostates.push_back(state);
               state.clear();

               qo.push_back(Quantum(2,0));

               //up down
               state = istate;
               state.push_crea_down(j);
               ostates.push_back(state);
               state.clear();

               qo.push_back(Quantum(1,1));

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

               qo.push_back(Quantum(1,1));

               //down down
               state = istate;
               state.push_crea_down(j);
               ostates.push_back(state);
               state.clear();

               qo.push_back(Quantum(0,2));

            }

            istate.clear();

         }

         //get merged row
         get_merged_index(-qi,q_merged,ind_merged,inverse);
         merged_row.push_back(inverse);

         //get merged column
         get_merged_index(qo,q_merged,ind_merged,inverse);
         merged_col.push_back(inverse);

         //first row
         for(int col = 0;col < ostates.size();++col){

            int ai = ostates[col].gsite(0);
            int as = ostates[col].gspin(0);

            int bi = ostates[col].gsite(1);
            int bs = ostates[col].gspin(1);

            int o = ij2o[no-1][no-1];
            int v = ab2v[ai-no][bi-no];

            int s = ov2s[o][v];

            if(as == 0 && bs == 1)
               push_anni_down_anni_up(s,no-1,0,col,1);
            else if(as == 1 && bs == 0)
               push_anni_down_anni_up(s,no-1,0,col,-1);

         }

         //singles coming in
         int row = 1;

         while(istates[row].size() == 1){

            for(int col = 0;col < ostates.size();++col)
               push_single_in_complement(no-1,istates[row],row,ostates[col],col);

            ++row;

         }

         //doubles coming in
         while(row < istates.size()){

            for(int col = 0;col < ostates.size();++col)
               push_double_complement(no-1,istates[row],row,ostates[col],col);

            ++row;

         }

      }
      else{//no < nv

         //last occupied: i = no - 1
         ostates.clear();
         qo.clear();

         qo.push_back(Quantum(1,1));//a_up a_down

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

         //get merged row
         get_merged_index(-qi,q_merged,ind_merged,inverse);
         merged_row.push_back(inverse);

         //get merged column
         get_merged_index(qo,q_merged,ind_merged,inverse);
         merged_col.push_back(inverse);

         istates = ostates;
         qi = qo;

         if(no + 1 == nv){

            ostates.clear();
            qo.clear();

            //first col closed
            state.push_id();
            ostates.push_back(state);
            state.clear();

            qo.push_back(Quantum::zero());

            //singles
            for(int i = no + 1;i < L;++i){

               state.push_crea_up(i);
               ostates.push_back(state);
               state.clear();

               qo.push_back(Quantum(1,0));

               state.push_crea_down(i);
               ostates.push_back(state);
               state.clear();

               qo.push_back(Quantum(0,1));

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

               qo.push_back(Quantum(1,1));

               for(int j = i + 1;j < L;++j){

                  //up up
                  state = istate;
                  state.push_crea_up(j);
                  ostates.push_back(state);
                  state.clear();

                  qo.push_back(Quantum(2,0));

                  //up down
                  state = istate;
                  state.push_crea_down(j);
                  ostates.push_back(state);
                  state.clear();

                  qo.push_back(Quantum(1,1));

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

                  qo.push_back(Quantum(1,1));

                  //down down
                  state = istate;
                  state.push_crea_down(j);
                  ostates.push_back(state);
                  state.clear();

                  qo.push_back(Quantum(0,2));

               }

               istate.clear();

            }

            //get merged row
            get_merged_index(-qi,q_merged,ind_merged,inverse);
            merged_row.push_back(inverse);

            //get merged column
            get_merged_index(qo,q_merged,ind_merged,inverse);
            merged_col.push_back(inverse);

            //first column closed: insert pairs
            for(int row = 0;row < istates.size();++row){

               //lets call in i,j
               int i = istates[row].gsite(1);
               int j = istates[row].gsite(0);

               int si = istates[row].gspin(1);
               int sj = istates[row].gspin(0);

               int o = ij2o[i][j];
               int v = ab2v[0][0];

               int s = ov2s[o][v];

               if(si == 0 && sj == 1)
                  push_crea_up_crea_down(s,no,row,0,1);
               else if(si == 1 && sj == 0)
                  push_crea_up_crea_down(s,no,row,0,-1);

            }

            //singles going out:
            int col = 1;

            while(ostates[col].size() == 1){

               for(int row = 0;row < istates.size();++row)
                  push_single_out_complement(no,istates[row],row,ostates[col],col);

               ++col;

            }

            //switch to outgoing pairs
            while(col < ostates.size()){

               for(int row = 0;row < istates.size();++row)
                  push_double_complement(no,istates[row],row,ostates[col],col);

               ++col;

            }

         }
         else{//nv > no + 1

            //for once the ostates and qo from previous site are correct for start of next site

            //all the virtuals: the have the signature of the operator they are going to: but opposite q-number
            for(int j = no + 1;j < L;++j){

               qo.push_back(Quantum(1,0));//go to crea up

               state.push_crea_up(j);
               ostates.push_back(state);
               state.clear();

               qo.push_back(Quantum(0,1));//go to crea down 

               state.push_crea_down(j);
               ostates.push_back(state);
               state.clear();

            }

            //last column complete closed:
            qo.push_back(Quantum(0,0));
            state.push_id();
            ostates.push_back(state);
            state.clear();

            //get merged row
            get_merged_index(-qi,q_merged,ind_merged,inverse);
            merged_row.push_back(inverse);

            //get merged column
            get_merged_index(qo,q_merged,ind_merged,inverse);
            merged_col.push_back(inverse);

            for(int col = istates.size();col < ostates.size() - 1;++col)
               for(int row = 0;row < istates.size();++row)
                  push_single_out_complement(no,istates[row],row,ostates[col],col);

            //last column: insert pairs
            for(int row = 0;row < istates.size();++row){

               //lets call in i,j
               int i = istates[row].gsite(1);
               int j = istates[row].gsite(0);

               int si = istates[row].gspin(1);
               int sj = istates[row].gspin(0);

               int o = ij2o[i][j];
               int v = ab2v[0][0];

               int s = ov2s[o][v];

               if(si == 0 && sj == 1)
                  push_crea_up_crea_down(s,no,row,ostates.size() - 1,1);
               else if(si == 1 && sj == 0)
                  push_crea_up_crea_down(s,no,row,ostates.size() - 1,-1);

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

                  qo.push_back(Quantum(1,0));//go to crea up

                  state.push_crea_up(j);
                  ostates.push_back(state);
                  state.clear();

                  qo.push_back(Quantum(0,1));//go to crea down 

                  state.push_crea_down(j);
                  ostates.push_back(state);
                  state.clear();

               }

               //last column closed
               state.push_id();
               ostates.push_back(state);
               state.clear();
               qo.push_back(Quantum(0,0));

               //get merged row
               get_merged_index(-qi,q_merged,ind_merged,inverse);
               merged_row.push_back(inverse);

               //get merged column
               get_merged_index(qo,q_merged,ind_merged,inverse);
               merged_col.push_back(inverse);

               int col = 0;

               //go past the pairs
               while(ostates[col].size() == 2)
                  ++col;

               //remaining virtuals
               while(col < ostates.size() - 1){

                  int row = 0;

                  while(istates[row].size() == 2){

                     push_single_out_complement(i,istates[row],row,ostates[col],col);

                     ++row;

                  }

                  ++col;

               }

               //last column: close down
               row = 0;

               while(istates[row].size() == 2){

                  //lets call in i,j
                  int ii = istates[row].gsite(1);
                  int ij = istates[row].gsite(0);

                  int si = istates[row].gspin(1);
                  int sj = istates[row].gspin(0);

                  int o = ij2o[ii][ij];
                  int v = ab2v[vind][vind];

                  int s = ov2s[o][v];

                  if(si == 0 && sj == 1)
                     push_crea_up_crea_down(s,i,row,ostates.size() - 1,1);
                  else if(si == 1 && sj == 0)
                     push_crea_up_crea_down(s,i,row,ostates.size() - 1,-1);

                  ++row;

               }

               istates = ostates;
               qi = qo;

            }

            ostates.clear();
            qo.clear();

            //first col closed
            state.push_id();
            ostates.push_back(state);
            state.clear();

            qo.push_back(Quantum::zero());

            //singles
            for(int i = nv + 1;i < L;++i){

               state.push_crea_up(i);
               ostates.push_back(state);
               state.clear();

               qo.push_back(Quantum(1,0));

               state.push_crea_down(i);
               ostates.push_back(state);
               state.clear();

               qo.push_back(Quantum(0,1));

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

               qo.push_back(Quantum(1,1));

               for(int j = i + 1;j < L;++j){

                  //up up
                  state = istate;
                  state.push_crea_up(j);
                  ostates.push_back(state);
                  state.clear();

                  qo.push_back(Quantum(2,0));

                  //up down
                  state = istate;
                  state.push_crea_down(j);
                  ostates.push_back(state);
                  state.clear();

                  qo.push_back(Quantum(1,1));

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

                  qo.push_back(Quantum(1,1));

                  //down down
                  state = istate;
                  state.push_crea_down(j);
                  ostates.push_back(state);
                  state.clear();

                  qo.push_back(Quantum(0,2));

               }

               istate.clear();

            }

            //get merged row
            get_merged_index(-qi,q_merged,ind_merged,inverse);
            merged_row.push_back(inverse);

            //get merged column
            get_merged_index(qo,q_merged,ind_merged,inverse);
            merged_col.push_back(inverse);

            //incoming pairs
            int row = 0;

            while(istates[row].size() == 2){

               //lets call in i,j
               int i = istates[row].gsite(1);
               int j = istates[row].gsite(0);

               int si = istates[row].gspin(1);
               int sj = istates[row].gspin(0);

               int o = ij2o[i][j];
               int v = ab2v[nv - no][nv - no];

               int s = ov2s[o][v];

               if(si == 0 && sj == 1)
                  push_crea_up_crea_down(s,nv,row,0,1);
               else if(si == 1 && sj == 0)
                  push_crea_up_crea_down(s,nv,row,0,-1);

               ++row;

            }

            //next the single columns
            int col = 1;

            while(ostates[col].size() == 1){

               row = 0;

               //doubles coming in
               while(istates[row].size() == 2){

                  push_single_out_complement(nv,istates[row],row,ostates[col],col);

                  ++row;

               }

               ++col;

            }

            while(col < ostates.size()){

               //transform from incoming to outgoing pairs
               row = 0;

               while(istates[row].size() == 2){

                  push_double_complement(nv,istates[row],row,ostates[col],col);

                  ++row;

               }

               ++col;

            }

         }

      }

   }

}

/**
 * copy constructor
 */
T2_2_mpo::T2_2_mpo(const T2_2_mpo &copy){ }

/**
 * destructor
 */
T2_2_mpo::~T2_2_mpo(){

   for(int i = 0;i < no;++i)
      delete [] ij2o[i];

   delete [] ij2o;

   for(int a = 0;a < nv;++a)
      delete [] ab2v[a];

   delete [] ab2v;

   for(int o = 0;o < no*(no + 1)/2;++o)
      delete [] ov2s[o];

   delete [] ov2s;

   delete [] list;

}

/**
 * add an anni down anni up to the list: on site with row and column and 
 */
void T2_2_mpo::push_anni_down_anni_up(int s,int site,int row,int col,int sign){

   std::vector<int> mpoind(6);

   mpoind[0] = site;
   mpoind[1] = row;//incoming virtual
   mpoind[2] = 0;//physical incoming
   mpoind[3] = 3;//physical outgoing
   mpoind[4] = col;//outgoing virtual
   mpoind[5] = sign;//sign

   list[s].push_back(mpoind);

}

/**
 * add an anni down anni up to the list: on site with row and column and 
 */
void T2_2_mpo::push_crea_up_crea_down(int s,int site,int row,int col,int sign){

   std::vector<int> mpoind(6);

   mpoind[0] = site;
   mpoind[1] = row;//incoming virtual
   mpoind[2] = 3;//physical incoming
   mpoind[3] = 0;//physical outgoing
   mpoind[4] = col;//outgoing virtual
   mpoind[5] = sign;//sign

   list[s].push_back(mpoind);

}

/**
 * get the complementary operator between in and out, when a single is coming in and a pair is going out and give the correct sign and operator to rearrange to normal order: for the T2 operator
 */
void T2_2_mpo::push_single_in_complement(int j,const Ostate &in,int row,const Ostate &out,int col){

   //lets call in i,j and out k,l
   int i = in.gsite(0);
   int si = in.gspin(0);

   int k = out.gsite(0);
   int sk = out.gspin(0);

   int l = out.gsite(1);
   int sl = out.gspin(1);

   //in virtual notation that is:
   int a = k - no;
   int b = l - no;

   int o,v,s;

   if(si == 0){//first site anni up spin

      if(sk == 0){//third site crea up spin

         if(sl == 0){//fourth crea up spin

            //+ij;ab
            o = ij2o[i][j];
            v = ab2v[a][b];
            s = ov2s[o][v];

            push_anni_up(s,j,row,col,1);

            //-ij;ba
            v = ab2v[b][a];
            s = ov2s[o][v];

            push_anni_up(s,j,row,col,-1);

         }
         else{//fourth crea down

            o = ij2o[i][j];
            v = ab2v[a][b];
            s = ov2s[o][v];

            push_anni_down_s(s,j,row,col,1);

         }

      }
      else{//a crea down

         if(sl == 0){//only one possible

            //-!
            o = ij2o[i][j];
            v = ab2v[b][a];
            s = ov2s[o][v];

            push_anni_down_s(s,j,row,col,-1);

         }

      }

   }
   else{//first site anni down

      if(sk == 0){//third site crea up

         if(sl == 1){//fourth site crea down

            //-!
            o = ij2o[i][j];
            v = ab2v[b][a];
            s = ov2s[o][v];

            push_anni_up(s,j,row,col,-1);

         }

      }
      else{//third site crea down

         if(sl == 1){//fourth site crea down

            //+ij;ab
            o = ij2o[i][j];
            v = ab2v[a][b];
            s = ov2s[o][v];

            push_anni_down_s(s,j,row,col,1);

            //-ij;ba
            v = ab2v[b][a];
            s = ov2s[o][v];

            push_anni_down_s(s,j,row,col,-1);

         }
         else{//fourth site crea up

            o = ij2o[i][j];
            v = ab2v[a][b];
            s = ov2s[o][v];

            push_anni_up(s,j,row,col,1);

         }

      }

   }

}

/**
 * add a crea up to the list: on site with row and column and 
 */
void T2_2_mpo::push_anni_up(int s,int site,int row,int col,int sign){

   std::vector<int> mpoind(6);

   mpoind[0] = site;
   mpoind[1] = row;//incoming virtual
   mpoind[4] = col;//outgoing virtual

   //0-> up
   mpoind[2] = 0;//physical incoming
   mpoind[3] = 1;//physical outgoing
   mpoind[5] = sign;//sign

   list[s].push_back(mpoind);

   //down-> up down
   mpoind[2] = 2;//physical incoming
   mpoind[3] = 3;//physical outgoing
   mpoind[5] = sign;//sign

   list[s].push_back(mpoind);

}

/**
 * add a crea up to the list: on site with row and column and 
 */
void T2_2_mpo::push_crea_up_s(int s,int site,int row,int col,int sign){

   std::vector<int> mpoind(6);

   mpoind[0] = site;
   mpoind[1] = row;//incoming virtual
   mpoind[4] = col;//outgoing virtual

   //0-> up
   mpoind[2] = 1;//physical incoming
   mpoind[3] = 0;//physical outgoing
   mpoind[5] = sign;//sign

   list[s].push_back(mpoind);

   //down-> up down
   mpoind[2] = 3;//physical incoming
   mpoind[3] = 2;//physical outgoing
   mpoind[5] = -sign;//sign

   list[s].push_back(mpoind);

}

/**
 * add a crea up to the list: on site with row and column and 
 */
void T2_2_mpo::push_anni_down_s(int s,int site,int row,int col,int sign){

   std::vector<int> mpoind(6);

   mpoind[0] = site;
   mpoind[1] = row;//incoming virtual
   mpoind[4] = col;//outgoing virtual

   //0-> down 
   mpoind[2] = 0;//physical incoming
   mpoind[3] = 2;//physical outgoing
   mpoind[5] = sign;//sign

   list[s].push_back(mpoind);

   //up-> up down
   mpoind[2] = 1;//physical incoming
   mpoind[3] = 3;//physical outgoing
   mpoind[5] = -sign;//sign

   list[s].push_back(mpoind);

}

/**
 * add a crea up to the list: on site with row and column and 
 */
void T2_2_mpo::push_crea_down(int s,int site,int row,int col,int sign){

   std::vector<int> mpoind(6);

   mpoind[0] = site;
   mpoind[1] = row;//incoming virtual
   mpoind[4] = col;//outgoing virtual

   //0-> down 
   mpoind[2] = 2;//physical incoming
   mpoind[3] = 0;//physical outgoing
   mpoind[5] = sign;//sign

   list[s].push_back(mpoind);

   //up-> up down
   mpoind[2] = 3;//physical incoming
   mpoind[3] = 1;//physical outgoing
   mpoind[5] = sign;//sign

   list[s].push_back(mpoind);

}

/**
 * get the value to attach to the operator between in and out
 * @return 0 if no connection, 1 if insert id and put value in teff
 */
void T2_2_mpo::push_double_complement(int site,const Ostate &in,int row,const Ostate &out,int col){

   //lets call in i,j and out l
   int i = in.gsite(1);
   int j = in.gsite(0);

   int si = in.gspin(1);
   int sj = in.gspin(0);

   int k = out.gsite(0);
   int sk = out.gspin(0);

   int l = out.gsite(1);
   int sl = out.gspin(1);

   //in virtual notation that is:
   int a = k - no;
   int b = l - no;

   int o,v,s;

   if(si == 0 && sj == 0){//anni up anni up coming in

      //only coupling with crea up crea up
      if(sk == 0 && sl == 0){

         //+ij;ab
         o = ij2o[i][j];
         v = ab2v[a][b];
         s = ov2s[o][v];

         push_id(s,site,row,col,1);

         //-ij;ba
         v = ab2v[b][a];
         s = ov2s[o][v];

         push_id(s,site,row,col,-1);

      }

   }
   else if(si == 0 && sj == 1){//anni up anni down coming in

      //coupling with up down and down up
      if(sk == 0 && sl == 1){

         //+ij;ab
         o = ij2o[i][j];
         v = ab2v[a][b];
         s = ov2s[o][v];

         push_id(s,site,row,col,1);

      }
      else if(sk == 1 && sl == 0){

         //-ij;ba
         o = ij2o[i][j];
         v = ab2v[b][a];
         s = ov2s[o][v];

         push_id(s,site,row,col,-1);

      }

   }
   else if(si == 1 && sj == 0){

      //coupling with up down and down up
      if(sk == 1 && sl == 0){

         //+ij;ab
         o = ij2o[i][j];
         v = ab2v[a][b];
         s = ov2s[o][v];

         push_id(s,site,row,col,1);

      }
      else if(sk == 0 && sl == 1){

         //-ij;ba
         o = ij2o[i][j];
         v = ab2v[b][a];
         s = ov2s[o][v];

         push_id(s,site,row,col,-1);

      }

   }
   else{//down down

      //only coupling with crea down crea down
      if(sk == 1 && sl == 1){

         //+ij;ab
         o = ij2o[i][j];
         v = ab2v[a][b];
         s = ov2s[o][v];

         push_id(s,site,row,col,1);

         //-ij;ba
         v = ab2v[b][a];
         s = ov2s[o][v];

         push_id(s,site,row,col,-1);

      }

   }

}

/**
 * add an id to the list: on site with row and column and index s(ia)
 */
void T2_2_mpo::push_id(int s,int site,int row,int col,int sign){

   std::vector<int> mpoind(6);

   mpoind[0] = site;
   mpoind[1] = row;//incoming virtual
   mpoind[4] = col;//outgoing virtual

   //0-> 0
   mpoind[2] = 0;//physical incoming
   mpoind[3] = 0;//physical outgoing
   mpoind[5] = sign;//sign

   list[s].push_back(mpoind);

   //1->1
   mpoind[2] = 1;//physical incoming
   mpoind[3] = 1;//physical outgoing
   mpoind[5] = sign;//sign

   list[s].push_back(mpoind);

   //2->2
   mpoind[2] = 2;//physical incoming
   mpoind[3] = 2;//physical outgoing
   mpoind[5] = sign;//sign

   list[s].push_back(mpoind);

   //3->3
   mpoind[2] = 3;//physical incoming
   mpoind[3] = 3;//physical outgoing
   mpoind[5] = sign;//sign

   list[s].push_back(mpoind);

}

ostream &operator<<(ostream &output,const T2_2_mpo &list_p){

   for(int s = 0;s < list_p.no*(list_p.no - 1)/2 * list_p.nv * list_p.nv;++s){

      int o = list_p.s2ov[s][0];
      int v = list_p.s2ov[s][1];

      int i = list_p.o2ij[o][0];
      int j = list_p.o2ij[o][1];

      int a = list_p.v2ab[v][0];
      int b = list_p.v2ab[v][1];

      output << std::endl;
      output << "element t(" << i << "," << j << ";" << a << "," << b << ")" << std::endl;
      output << std::endl;

      std::vector<int> br;
      std::vector<int> bc;

      for(int ind = 0;ind < list_p.list[s].size();++ind){

         int site = list_p.list[s][ind][0];
         int r = list_p.list[s][ind][1];
         int pi = list_p.list[s][ind][2];
         int po = list_p.list[s][ind][3];
         int c = list_p.list[s][ind][4];
         int sign = list_p.list[s][ind][5];

         br = list_p.merged_row[site][r];
         bc = list_p.merged_col[site][c];

         output << site << "\tunmerged index: (" << r << "," << pi << "," << po << "," << c << ")\t" << 

            "merged block:\t(" << br[0] << "," << pi << "," << po << "," << bc[0] << ")\t" << 

            "merged block index:\t(" << br[1] << "," << 0 << "," << 0 << "," << bc[1] << ")" << 

            "\tsign = " << sign << std::endl;

      } 
   }

   return output;

}

/**
 * get the sum of the elements from input MPO corresponding to the t2 element t(i,j,a,b).
 */
template <class Q>
double T2_2_mpo::get(const MPO<Q> &grad,int i,int j,int a,int b,bool merged) const{

   int o = ij2o[i][j];
   int vi = ab2v[a][b];

   int s = ov2s[o][vi];

   std::vector<int> v;

   IVector<4> ind;
   QSDArray<4>::const_iterator it;

   if(merged){

      std::vector<int> br;
      std::vector<int> bc;

      double sum = 0.0;

      for(int n = 0;n < list[s].size();++n){

         v = list[s][n];

         br = merged_row[v[0]][v[1]];
         bc = merged_col[v[0]][v[4]];

         ind[0] = br[0];
         ind[1] = v[2];
         ind[2] = v[3];
         ind[3] = bc[0];

         it = grad[v[0]].find(ind);

         if(it != grad[v[0]].end())//if the block exists!
            sum += (*it->second).at(br[1],0,0,bc[1]) * v[5];

      }

      return sum;

   }
   else{

      double sum = 0.0;

      for(int n = 0;n < list[s].size();++n){

         v = list[s][n];

         ind[0] = v[1];
         ind[1] = v[2];
         ind[2] = v[3];
         ind[3] = v[4];

         it = grad[v[0]].find(ind);

         if(it != grad[v[0]].end())//if the block exists!
            sum += (*it->second).at(0,0,0,0) * v[5];

      }

      return sum;

   }

}

/**
 * print the elements from input MPO corresponding to the t2 element t(i,j,a,b), for testing
 */
template <class Q>
void T2_2_mpo::print(const MPO<Q> &grad,int i,int j,int a,int b,bool merged) const {

   int o = ij2o[i][j];
   int vi = ab2v[a][b];

   int s = ov2s[o][vi];

   std::vector<int> v;

   IVector<4> ind;
   QSDArray<4>::const_iterator it;

   if(merged){

      std::vector<int> br;
      std::vector<int> bc;

      double sum = 0.0;

      for(int n = 0;n < list[s].size();++n){

         v = list[s][n];

         br = merged_row[v[0]][v[1]];
         bc = merged_col[v[0]][v[4]];

         ind[0] = br[0];
         ind[1] = v[2];
         ind[2] = v[3];
         ind[3] = bc[0];

         it = grad[v[0]].find(ind);

         if(it != grad[v[0]].end())//if the block exists!
            cout << n << "\t" << (*it->second).at(br[1],0,0,bc[1]) * v[5] << endl;

      }

   }
   else{

      double sum = 0.0;

      for(int n = 0;n < list[s].size();++n){

         v = list[s][n];

         ind[0] = v[1];
         ind[1] = v[2];
         ind[2] = v[3];
         ind[3] = v[4];

         it = grad[v[0]].find(ind);

         if(it != grad[v[0]].end())//if the block exists!
            cout << n << "\t" << (*it->second).at(0,0,0,0) * v[5] << endl;

      }

   }

}

/**
 * get the complementary operator between in and out, when a pair is coming and one is going out and give the correct sign to rearrange to normal order: for the T2 operator
 */
void T2_2_mpo::push_single_out_complement(int site,const Ostate &in,int row,const Ostate &out,int col){

   //lets call in i,j and out l
   int i = in.gsite(1);
   int j = in.gsite(0);

   int si = in.gspin(1);
   int sj = in.gspin(0);

   int l = out.gsite(0);
   int sl = out.gspin(0);

   //in virtual notation that is:
   int a = site - no;
   int b = l - no;

   int o,v,s;

   if(si == 0){//first site anni up spin

      if(sj == 0){//second site anni up spin

         if(sl == 0){//only coupling with crea up spin

            //+ij;ab
            o = ij2o[i][j];
            v = ab2v[a][b];
            s = ov2s[o][v];

            push_crea_up_s(s,site,row,col,1);

            //-ij;ba
            v = ab2v[b][a];
            s = ov2s[o][v];

            push_crea_up_s(s,site,row,col,-1);

         }

      }
      else{//second site anni down spin

         if(sl == 0){//out crea up

            //-ij;ba
            o = ij2o[i][j];
            v = ab2v[b][a];
            s = ov2s[o][v];

            push_crea_down(s,site,row,col,-1);

         }
         else{//out crea down

            //ij;ab
            o = ij2o[i][j];
            v = ab2v[a][b];
            s = ov2s[o][v];

            push_crea_up_s(s,site,row,col,1);

         }

      }

   }
   else{//first site anni down

      if(sj == 0){//second site anni up

         if(sl == 0){//out crea up

            //ij;ab
            o = ij2o[i][j];
            v = ab2v[a][b];
            s = ov2s[o][v];

            push_crea_down(s,site,row,col,1);

         }
         else{//out crea down

            o = ij2o[i][j];
            v = ab2v[b][a];
            s = ov2s[o][v];

            push_crea_up_s(s,site,row,col,-1);

         }

      }
      else{//second site anni down

         //only coupling with crea down
         if(sl == 1){

            //+ij;ab
            o = ij2o[i][j];
            v = ab2v[a][b];
            s = ov2s[o][v];

            push_crea_down(s,site,row,col,1);

            //-ij;ba
            v = ab2v[b][a];
            s = ov2s[o][v];

            push_crea_down(s,site,row,col,-1);

         }

      }

   }

}

template double T2_2_mpo::get(const MPO<Quantum> &,int,int,int,int,bool) const;
template void T2_2_mpo::print(const MPO<Quantum> &,int,int,int,int,bool) const;
