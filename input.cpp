#include <iostream>
#include <fstream>
#include <iomanip>

using std::cout;
using std::endl;
using std::ostream;
using std::ifstream;

#include "include.h"

/**
 * read in the one electron integrals from a file, to a Darray<2> object, giving a list of the orbitals orderings
 */
void read_oei(const char *filename,DArray<2> &t,const std::vector<int> &order){

   ifstream input(filename);

   int i,j;
   double value;

   while(input >> i >> j >> value){

      t(order[i],order[j]) = value;
      t(order[j],order[i]) = value;

   }

}

/**
 * read in the two electron integrals from a file, to a Darray<2> object, giving a list of the orbitals orderings
 */
void read_tei(const char *filename,DArray<4> &V,const std::vector<int> &order){

   ifstream input(filename);

   int a,b,c,d;
   double value;

   while(input >> a >> c >> b >> d >> value){

      V(order[a],order[b],order[c],order[d]) = value;//ab;cd
      V(order[b],order[a],order[d],order[c]) = value;//ba;dc
      V(order[c],order[b],order[a],order[d]) = value;//cb;ad
      V(order[b],order[c],order[d],order[a]) = value;//bc;da
      V(order[a],order[d],order[c],order[b]) = value;//ad;cb
      V(order[d],order[a],order[b],order[c]) = value;//da;bc
      V(order[c],order[d],order[a],order[b]) = value;//cd;ab
      V(order[d],order[c],order[b],order[a]) = value;//dc;ba

   }

}
