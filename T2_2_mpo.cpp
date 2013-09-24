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
   vector<int> v(2);

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
      for(int vi = 0;vi < no*(no + 1)/2;++vi){

         ov2s[o][vi] = s;
         ++s;

         v[0] = 0;
         v[1] = vi;

         s2ov.push_back(v);

      }

   list = new std::vector< vector<int> > [s2ov.size()];

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

   for(int i = 1;i < no - 1;++i){

      ostates.clear();
      ostates.push_back(istates[0]);

      state.push_anni_down(i);
      ostates.push_back(state);
      state.clear();

      state.push_anni_up(i);
      ostates.push_back(state);
      state.clear();

      int row = 1;

      while(istates[row].size() == 1){

         ostates.push_back(istates[row]);
         ++row;

      }

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

         ++row;

      }

      //add an up
      row = 1;

      while(istates[row].size() == 1){

         state.push_anni_up(i);
         state.insert(state.end(),istates[row].begin(),istates[row].end());
         ostates.push_back(state);
         state.clear();

         ++row;

      }

      //id for the pairs coming in
      while(row < istates.size()){

         ostates.push_back(istates[row]);
         ++row;

      }

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
