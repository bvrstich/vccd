#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

using std::endl;
using std::ostream;
using std::ofstream;
using std::ifstream;

#include "include.h"

std::vector< QSDArray<2> ** > Operator::op;
std::vector< bool * > Operator::sparse;
std::vector< std::vector<int> > Operator::dim;

int Operator::L;

//!standard constructor
Operator::Operator(){ }

//!copy constructor
Operator::Operator(const Operator &){ }

//!virtual destructor
Operator::~Operator(){ } 

/**
 * initialize all the operators, construct the complementary operators for the qc Hamiltonian
 */
void Operator::init(const DArray<2> &t,const DArray<4> &V){ 

   L = t.shape(0);

   dim.resize(L);

   dim[0].resize(2);

   //row dim
   dim[0][0] = 1;

   for(int i = 1;i < L;++i){

      dim[i - 1][1] = HamOp::ostates[i - 1].size();

      dim[i].resize(2);

      dim[i][0] = dim[i - 1][1];

   }

   dim[L - 1][1] = 1;

   sparse.resize(L);
   op.resize(L);

   for(int i = 0;i < L;++i){

      sparse[i] = new bool [dim[i][0]*dim[i][1]];
      op[i] = new QSDArray<2> * [dim[i][0]*dim[i][1]];

   }

      /*
         Qshapes<Quantum> qp;
         physical(qp);

         Quantum qt;
         DArray<2> Ip(1,1);
         DArray<2> Im(1,1);

         Ip = 1.0;
         Im = -1.0;

      //first id
      qt = Quantum::zero();

      id.resize(qt,make_array(qp,-qp));

      id.insert(shape(0,0),Ip);
      id.insert(shape(1,1),Ip);
      id.insert(shape(2,2),Ip);
      id.insert(shape(3,3),Ip);

      //sign
      s.resize(qt,make_array(qp,-qp));

      Ip = 1.0;

      s.insert(shape(0,0),Ip);
      s.insert(shape(1,1),Im);
      s.insert(shape(2,2),Im);
      s.insert(shape(3,3),Ip);

      //create up
      qt = Quantum(1,0);

      cu.resize(qt,make_array(qp,-qp));

      cu.insert(shape(1,0),Ip);
      cu.insert(shape(3,2),Ip);

      //create up sign
      cus.resize(qt,make_array(qp,-qp));

      cus.insert(shape(1,0),Ip);
      cus.insert(shape(3,2),Im);

      //create down
      qt = Quantum(0,1);

      cd.resize(qt,make_array(qp,-qp));

      cd.insert(shape(2,0),Ip);
      cd.insert(shape(3,1),Ip);

      //create down sign
      cds.resize(qt,make_array(qp,-qp));

      cds.insert(shape(2,0),Ip);
      cds.insert(shape(3,1),Im);

      //annihilator of up spin
      qt = Quantum(-1,0);

      au.resize(qt,make_array(qp,-qp));

      au.insert(shape(0,1),Ip);
      au.insert(shape(2,3),Ip);

      //annihilator of up spin sign
      aus.resize(qt,make_array(qp,-qp));

      aus.insert(shape(0,1),Ip);
      aus.insert(shape(2,3),Im);

   //annihilator of down spin
   qt = Quantum(0,-1);

   ad.resize(qt,make_array(qp,-qp));

   ad.insert(shape(0,2),Ip);
   ad.insert(shape(1,3),Ip);

   //annihilator of down spin with sign
   ads.resize(qt,make_array(qp,-qp));

   ads.insert(shape(0,2),Ip);
   ads.insert(shape(1,3),Im);

   //creator of an up down pair on site
   qt = Quantum(1,1);

   cucd.resize(qt,make_array(qp,-qp));

   cucd.insert(shape(3,0),Ip);

   //create up annihilate up
   qt = Quantum::zero();

   cuau.resize(qt,make_array(qp,-qp));

   cuau.insert(shape(1,1),Ip);
   cuau.insert(shape(3,3),Ip);

   //create up annihilate down
   qt = Quantum(1,-1);

   cuad.resize(qt,make_array(qp,-qp));

   cuad.insert(shape(1,2),Ip);

   //create down annihilate up
   qt = Quantum(-1,1);

   cdau.resize(qt,make_array(qp,-qp));

   cdau.insert(shape(2,1),Ip);

   //create down annihilate down --> n_down
   qt = Quantum::zero();

   cdad.resize(qt,make_array(qp,-qp));

   cdad.insert(shape(2,2),Ip);
   cdad.insert(shape(3,3),Ip);

   //insert annihilate up annihilate down: extra minus sign! because it should be annihilate down annihilate up
   qt = Quantum(-1,-1);

   auad.resize(qt,make_array(qp,-qp));

   auad.insert(shape(0,3),Im);

   //create the complementary triple operators
   tcuf.resize(L - 1);
   tcdf.resize(L - 1);
   tadf.resize(L - 1);
   tauf.resize(L - 1);

   for(int i = 0;i < L-1;++i){

      tcuf[i].resize(L - 1 - i);
      tcdf[i].resize(L - 1 - i);
      tadf[i].resize(L - 1 - i);
      tauf[i].resize(L - 1 - i);

      for(int j = i + 1;j < L;++j){

         //complementary of create up, first: triple comes to the left
         qt = Quantum(-1,0);

         tcuf[i][j - i - 1].resize(qt,make_array(qp,-qp));

         Ip = t(i,j);
         tcuf[i][j - i - 1].insert(shape(0,1),Ip);

         Ip = -t(i,j) - V(i,i,j,i);
         tcuf[i][j - i - 1].insert(shape(2,3),Ip);

         //complementary of create down, first: triple comes to the left
         qt = Quantum(0,-1);

         tcdf[i][j - i - 1].resize(qt,make_array(qp,-qp));

         Ip = t(i,j);
         tcdf[i][j - i - 1].insert(shape(0,2),Ip);

         Ip = t(i,j) + V(i,i,j,i);
         tcdf[i][j - i - 1].insert(shape(1,3),Ip);

         //complementary of anni down, first: triple comes to the left
         qt = Quantum(0,1);

         tadf[i][j - i - 1].resize(qt,make_array(qp,-qp));

         Ip = t(i,j);
         tadf[i][j - i - 1].insert(shape(2,0),Ip);

         Ip = t(i,j) + V(i,i,j,i);
         tadf[i][j - i - 1].insert(shape(3,1),Ip);

         //complementary of anni up, first: triple comes to the left
         qt = Quantum(1,0);

         tauf[i][j - i - 1].resize(qt,make_array(qp,-qp));

         Ip = t(i,j);
         tauf[i][j - i - 1].insert(shape(1,0),Ip);

         Ip = -t(i,j) - V(i,i,j,i);
         tauf[i][j - i - 1].insert(shape(2,3),Ip);

      }

   }

   local.resize(L);

   for(int i = 0;i < L;++i){

      //local term
      qt = Quantum::zero();

      local[i].resize(qt,make_array(qp,-qp));

      Ip = t(i,i);

      local[i].insert(shape(1,1),Ip);
      local[i].insert(shape(2,2),Ip);

      Ip = 2*t(i,i) + V(i,i,i,i);

      local[i].insert(shape(3,3),Ip);

   }
   */
}

/**
 * delete all the lists and stuff
 */
void Operator::clear(){

   for(int i = 0;i < L;++i)
      delete [] sparse[i];

   for(int i = 0;i < L;++i)
      delete [] op[i];

}
