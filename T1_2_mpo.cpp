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

   int L = no + nv;

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

      while(row < istates.size()){

         ostates.push_back(istates[row]);

         ++row;

      }

      istates = ostates;

   }

   //final occupied
   ostates.clear();

   state.push_anni_down(no - 1);
   ostates.push_back(state);
   state.clear();

   state.push_anni_up(no - 1);
   ostates.push_back(state);
   state.clear();

   int row = 1;

   while(row < istates.size()){

      ostates.push_back(istates[row]);

      ++row;

   }

   istates = ostates;

   if(no == nv){

      ostates.clear();

      state.push_id();
      ostates.push_back(state);
      state.clear();

      for(int i = no + 1;i < L;++i){

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

         if(si == 0)
            push_crea_up(ia2s[ii][0],no,row,0);
         else
            push_crea_down_s(ia2s[ii][0],no,row,0);

      }

      //insert the rest of the matrix: t's!
      for(int row = 0;row < istates.size();++row){

         int ii = istates[row].gsite(0);
         int si = istates[row].gspin(0);

         for(int col = 1;col < ostates.size();++col){

            int ia = ostates[col].gsite(0);
            int sa = ostates[col].gspin(0);

            if(si == sa)
               push_sign(ia2s[ii][ia - no],no,row,col);

         }

      }

      istates = ostates;

   }
   else{//no < nv

      //first virtual: first no states are the same: so no change in ostates and qo, only id added at the end for closing terms
      state.push_id();
      ostates.push_back(state);
      state.clear();

      //insert the t's for the last column
      row = 0;
      int column = istates.size();

      while(row < istates.size()){

         int j = istates[row].gsite(0);
         int sj = istates[row].gspin(0);

         if(sj == 0)
            push_crea_up(ia2s[j][0],no,row,column);
         else
            push_crea_down_s(ia2s[j][0],no,row,column);

         ++row;

      }

      istates = ostates;

      //all virtuals until nv
      for(int i = no + 1;i < nv - 1;++i){

         int vind = i - no;

         //insert the t's for the last column
         row = 0;
         column = istates.size() - 1;

         while(row < istates.size() - 1){

            int j = istates[row].gsite(0);
            int sj = istates[row].gspin(0);

            //fermion sign!
            if(sj == 0)
               push_crea_up(ia2s[j][vind],i,row,column);
            else
               push_crea_down_s(ia2s[j][vind],i,row,column);

            ++row;

         }

      }

      istates = ostates;

      ostates.clear();

      //switch to virtuals outgoing
      int vind = nv - no - 1;

      state.push_id();
      ostates.push_back(state);
      state.clear();

      for(int i = nv;i < L;++i){

         state.push_crea_down(i);
         ostates.push_back(state);
         state.clear();

         state.push_crea_up(i);
         ostates.push_back(state);
         state.clear();

      }

      //first column is closed
      for(int row = 0;row < istates.size()-1;++row){

         int ii = istates[row].gsite(0);
         int si = istates[row].gspin(0);

         if(si == 0)
            push_crea_up(ia2s[ii][vind],nv-1,row,0);
         else
            push_crea_down_s(ia2s[ii][vind],nv-1,row,0);

      }

      //insert the rest of the matrix: t's!
      for(int row = 0;row < istates.size();++row){

         int ii = istates[row].gsite(0);
         int si = istates[row].gspin(0);

         for(int col = 1;col < ostates.size();++col){

            int ia = ostates[col].gsite(0);
            int sa = ostates[col].gspin(0);

            if(si == sa)
               push_sign(ia2s[ii][ia - no],nv-1,row,col);

         }

      }

   }

}

/**
 * copy constructor
 */
T1_2_mpo::T1_2_mpo(const T1_2_mpo &copy){ }

/**
 * destructor
 */
T1_2_mpo::~T1_2_mpo(){ 

   for(int i = 0;i < no;++i)
      delete [] ia2s[i];

   delete [] ia2s;

   for(int s = 0;s < no*nv;++s)
      delete [] s2ia[s];

   delete [] s2ia;

   delete [] list;
   
}

/**
 * add a crea up to the list: on site with row and column and 
 */
void T1_2_mpo::push_crea_up(int s,int site,int row,int col){

   std::vector<int> mpoind(6);

   mpoind[0] = site;
   mpoind[1] = row;//incoming virtual
   mpoind[4] = col;//outgoing virtual

   //0-> up
   mpoind[2] = 1;//physical incoming
   mpoind[3] = 0;//physical outgoing
   mpoind[5] = 1;//sign

   list[s].push_back(mpoind);

   //down-> up down
   mpoind[2] = 3;//physical incoming
   mpoind[3] = 2;//physical outgoing
   mpoind[5] = 1;//sign

   list[s].push_back(mpoind);

}

/**
 * add a crea up to the list: on site with row and column and 
 */
void T1_2_mpo::push_crea_down_s(int s,int site,int row,int col){

   std::vector<int> mpoind(6);

   mpoind[0] = site;
   mpoind[1] = row;//incoming virtual
   mpoind[4] = col;//outgoing virtual

   //0-> down 
   mpoind[2] = 2;//physical incoming
   mpoind[3] = 0;//physical outgoing
   mpoind[5] = 1;//sign

   list[s].push_back(mpoind);

   //up-> up down
   mpoind[2] = 3;//physical incoming
   mpoind[3] = 1;//physical outgoing
   mpoind[5] = -1;//sign

   list[s].push_back(mpoind);

}

/**
 * add an id to the list: on site with row and column and index s(ia)
 */
void T1_2_mpo::push_sign(int s,int site,int row,int col){

   std::vector<int> mpoind(6);

   mpoind[0] = site;
   mpoind[1] = row;//incoming virtual
   mpoind[4] = col;//outgoing virtual

   //0-> 0
   mpoind[2] = 0;//physical incoming
   mpoind[3] = 0;//physical outgoing
   mpoind[5] = 1;//sign

   list[s].push_back(mpoind);

   //1->1
   mpoind[2] = 1;//physical incoming
   mpoind[3] = 1;//physical outgoing
   mpoind[5] = -1;//sign

   list[s].push_back(mpoind);

   //2->2
   mpoind[2] = 2;//physical incoming
   mpoind[3] = 2;//physical outgoing
   mpoind[5] = -1;//sign

   list[s].push_back(mpoind);

   //3->3
   mpoind[2] = 3;//physical incoming
   mpoind[3] = 3;//physical outgoing
   mpoind[5] = 1;//sign

   list[s].push_back(mpoind);

}


ostream &operator<<(ostream &output,const T1_2_mpo &list_p){

   for(int s = 0;s < list_p.no*list_p.nv;++s){

      int i = list_p.s2ia[s][0];
      int a = list_p.s2ia[s][1];
      
      output << std::endl;
      output << "element t(" << i << "," << a << ")" << std::endl;
      output << std::endl;

      for(int ind = 0;ind < list_p.list[s].size();++ind){

         output << list_p.list[s][ind][0] << "\t(" << list_p.list[s][ind][1] << "," << list_p.list[s][ind][2] << "," << list_p.list[s][ind][3] << 
         
            "," << list_p.list[s][ind][4] << ")\tsign = " << list_p.list[s][ind][5] << std::endl;

      }

   }

   return output;

}

/**
 * get the sum of the elements from input MPO corresponding to the t1 element t(i,a).
 */
template <class Q>
double T1_2_mpo::get(const MPO<Q> &grad,int i,int a){

   int s = ia2s[i][a];

   vector<int> v;
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

template double T1_2_mpo::get(const MPO<Quantum> &,int,int);
