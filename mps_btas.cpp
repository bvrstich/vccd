#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstdlib>

using std::cout;
using std::endl;
using std::ofstream;
using std::ifstream;

#include "FermiQuantum.h"
namespace btas { typedef FermiQuantum Quantum; }; // Defined as default quantum number class

/**
 * simple random number generator
 */
double rgen() { return 2.0*(static_cast<double>(rand())/RAND_MAX) - 1.0; }

#include "include.h"

using namespace btas;
using namespace mps;

int main(void){

   cout.precision(10);
   srand(time(NULL));

   //lenght of the chain
   int L = 14;

   //number of particles
   int n_u = 1;
   int n_d = 1;

   int no = n_u;
   int nv = L - no;

   Ostate::construct_oplist(L);

   Qshapes<Quantum> qp;
   physical(qp);

   MPS<Quantum> X = create(L,Quantum(n_u,n_d),qp,100,rgen);
   normalize(X);
   MPS<Quantum> Y = create(L,Quantum(n_u,n_d),qp,100,rgen);
   normalize(Y);

   DArray<2> t(L,L);

   for(int i = 0;i < L;++i)
      for(int j = 0;j < L;++j){

         double value = rgen();

         t(i,j) = value;
         t(j,i) = value;

      }


   MPO<Quantum> T = one_body<Quantum>(t);
   compress(T,mps::Right,0);
   compress(T,mps::Left,0);

   MPS<Quantum> TX = T*X;

   SDArray<1> S;//singular values
   QSDArray<2> V;//V^T
   QSDArray<3> U;//U --> unitary left normalized matrix

   for(int i = 0;i < L - 1;++i){

      //then svd
      QSDgesvd(RightArrow,TX[i],S,U,V,0);

      //copy unitary to mpx
      QSDcopy(U,TX[i]);

      //paste S and V together
      SDdidm(S,V);

      //and multiply with mpx on the next site
      U = TX[i + 1];

      //when compressing dimensions will change, so reset:
      TX[i + 1].clear();

      QSDcontract(1.0,V,shape(1),U,shape(0),0.0,TX[i + 1]);

      cout << i << "\t" << dot(mps::Left,TX,TX) << endl;

   }

   int i = L - 1;

   cout << endl;
   cout << TX[i].qshape() << endl;
   cout << TX[i].dshape() << endl;
   cout << endl;

   //then svd
   QSDgesvd(RightArrow,TX[i],S,U,V,0);

   //reverse the svd
   SDdimd(U,S);

   QSDArray<3> tmp;

   QSDcontract(1.0,U,shape(2),V,shape(0),0.0,tmp);

   QSDaxpy(-1.0,TX[i],tmp);
   cout << QSDdotc(tmp,tmp) << endl;

/*
   //copy unitary to mpx
   QSDcopy(U,TX[i]);

   //paste S and V together
   SDdidm(S,V);

   //and multiply with mpx on the next site
   U = TX[i + 1];

   //when compressing dimensions will change, so reset:
   TX[i + 1].clear();

   QSDcontract(1.0,V,shape(1),U,shape(0),0.0,TX[i + 1]);

   cout << i << "\t" << dot(mps::Left,TX,TX) << endl;
*/
   /*
   //make the HF state
   std::vector<int> occ(L);

   for(int i = 0;i < no;++i)
   occ[i] = 3;

   for(int i = no;i < L;++i)
   occ[i] = 0;

   MPS<Quantum> hf = product_state(L,qp,occ);
   MPS<Quantum> HA = qc*hf;
    */
   return 0;

}
