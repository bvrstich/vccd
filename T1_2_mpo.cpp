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

         ++s

      }

   list = new std::vector<int*> [no*nv];


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
