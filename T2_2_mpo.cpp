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

      istates = ostates;

   }

   if(no == nv){

      //last occupied: i = no - 1: switch from incoming to outgoing
      ostates.clear();

      Ostate istate;

      for(int i = no;i < L;++i){
         
         //first up
         istate.push_crea_up(i);

         state = istate;
         state.push_crea_down(i);
         ostates.push_back(state);
         state.clear();

         for(int j = i + 1;j < L;++j){

            //up up
            state = istate;
            state.push_crea_up(j);
            ostates.push_back(state);
            state.clear();

            //up down
            state = istate;
            state.push_crea_down(j);
            ostates.push_back(state);
            state.clear();

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

            //down down
            state = istate;
            state.push_crea_down(j);
            ostates.push_back(state);
            state.clear();

         }

         istate.clear();

      }

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
            push_single_complement(no-1,istates[row],row,ostates[col],col);

         ++row;

      }

      //doubles coming in
      while(row < istates.size()){

         for(int col = 0;col < ostates.size();++col)
            push_double_complement(no-1,istates[row],row,ostates[col],row);

         ++row;

      }

      istates = ostates;

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
 * get the complementary operator between in and out, when a single is coming in and a pair is going out and give the correct sign and operator to rearrange to normal order: for the T2 operator
 */
void T2_2_mpo::push_single_complement(int j,const Ostate &in,int row,const Ostate &out,int col){

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

      for(int ind = 0;ind < list_p.list[s].size();++ind){

         output << list_p.list[s][ind][0] << "\t(" << list_p.list[s][ind][1] << "," << list_p.list[s][ind][2] << "," << list_p.list[s][ind][3] << 
         
            "," << list_p.list[s][ind][4] << ")\tsign = " << list_p.list[s][ind][5] << std::endl;

      }

   }

   return output;

}

/**
 * get the sum of the elements from input MPO corresponding to the t2 element t(i,j,a,b).
 */
template <class Q>
double T2_2_mpo::get(const MPO<Q> &grad,int i,int j,int a,int b){

   int o = ij2o[i][j];
   int vi = ab2v[a][b];

   int s = ov2s[o][vi];

   std::vector<int> v;

   IVector<4> ind;
   QSDArray<4>::const_iterator it;

   double sum = 0.0;

   for(int n = 0;n < list[s].size();++n){

      v = list[s][n];

      ind[0] = v[1];
      ind[1] = v[2];
      ind[2] = v[3];
      ind[3] = v[4];

      it = grad[v[0]].find(ind);
      sum += (*it->second).at(0,0,0,0) * v[5];

   }

   return sum;

}

template double T2_2_mpo::get(const MPO<Quantum> &,int,int,int,int);
