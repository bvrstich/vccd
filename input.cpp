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

   t = 0.0;

   while(input >> value >> i >> j){

      t(order[i - 1],order[j - 1]) = value;
      t(order[j - 1],order[i - 1]) = value;

   }

}

/**
 * Fill the one electron integral matrix random
 */
void random_oei(DArray<2> &t){

   int L = t.shape(0);

   t = 0.0;

   for(int i = 0;i < L;++i)
      for(int j = i;j < L;++j){

         double value = 2.0*(static_cast<double>(rand())/RAND_MAX) - 1.0;

         t(i,j) = value;
         t(j,i) = value;

      }

}

/**
 * read in the two electron integrals from a file, to a Darray<2> object, giving a list of the orbitals orderings
 */
void read_tei(const char *filename,DArray<4> &V,const std::vector<int> &order){

   ifstream input(filename);

   int a,b,c,d;
   double value;

   V = 0.0;

   while(input >> value >> a >> c >> b >> d){

      V(order[a-1],order[b-1],order[c-1],order[d-1]) = value;
      V(order[b-1],order[a-1],order[d-1],order[c-1]) = value;
      V(order[c-1],order[b-1],order[a-1],order[d-1]) = value;
      V(order[b-1],order[c-1],order[d-1],order[a-1]) = value;
      V(order[a-1],order[d-1],order[c-1],order[b-1]) = value;
      V(order[d-1],order[a-1],order[b-1],order[c-1]) = value;
      V(order[c-1],order[d-1],order[a-1],order[b-1]) = value;
      V(order[d-1],order[c-1],order[b-1],order[a-1]) = value;

   }

}

/**
 * Fill the two-electron matrix whit random number
 */
void random_tei(DArray<4> &V){

   int L = V.shape(0);

   V = 0.0;

   for(int a = 0;a < L;++a)
      for(int b = 0;b < L;++b)
         for(int c = 0;c < L;++c)
            for(int d = 0;d < L;++d){

               double value = 2.0*(static_cast<double>(rand())/RAND_MAX) - 1.0;

               V(a,b,c,d) = value;
               V(b,a,d,c) = value;
               V(c,b,a,d) = value;
               V(b,c,d,a) = value;
               V(a,d,c,b) = value;
               V(d,a,b,c) = value;
               V(c,d,a,b) = value;
               V(d,c,b,a) = value;

            }


}
