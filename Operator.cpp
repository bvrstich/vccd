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
Qshapes<Quantum> Operator::qp;

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
   
   qp.clear();

   qp.push_back(Quantum(0,0));//zero
   qp.push_back(Quantum(1,0));//up
   qp.push_back(Quantum(0,1));//down
   qp.push_back(Quantum(1,1));//up down pair

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

      //initialize sparse on all false
      for(int row = 0;row < dim[i][0];++row)
         for(int col = 0;col < dim[i][1];++col)
            sparse[i][row*dim[i][1] + col] = false;

   }

   //now make the necessary operators!
   int row = 0;
   int col = 0;

   //construct id
   construct_id(0,0,0,1.0);

   //construct singles
   construct_cus(0,0,1,1.0);
   construct_cd(0,0,2,1.0);
   construct_aus(0,0,3,1.0);
   construct_ad(0,0,4,1.0);

   col = 5;

   //construct doubles
   construct_cucd(0,0,5,1.0);
   construct_cuau(0,0,6,1.0);
   construct_cuad(0,0,7,1.0);
   construct_cdau(0,0,8,1.0);
   construct_cdad(0,0,9,1.0);
   construct_adau(0,0,10,1.0);

   col = 11;

   //construct triples: complementary operator
   for(int j = 1;j < L;++j){

      construct_tcuf(0,0,col,t(0,j),V(0,0,j,0));col++;
      construct_tcdf(0,0,col,t(0,j),V(0,0,j,0));col++;
      construct_tauf(0,0,col,t(0,j),V(0,0,j,0));col++;
      construct_tadf(0,0,col,t(0,j),V(0,0,j,0));col++;

   }

   //last term:
   construct_local(0,0,col,t(0,0),V(0,0,0,0));

   //middle tensors
   for(int i = 1;i < L/2;++i){

      row = 0;
      col = 0;

      //construct id
      construct_id(i,0,0,1.0);

      //construct singles
      construct_cus(i,0,1,1.0);
      construct_cd(i,0,2,1.0);
      construct_aus(i,0,3,1.0);
      construct_ad(i,0,4,1.0);

      //construct signs
      row = 1;
      col = 5;

      while(HamOp::ostates[i-1][row].size() == 1){

         construct_s(i,row,col,1.0);
         ++row;
         ++col;

      }

      //construct doubles
      construct_cucd(i,0,col,1.0);++col;
      construct_cuau(i,0,col,1.0);++col;
      construct_cuad(i,0,col,1.0);++col;
      construct_cdau(i,0,col,1.0);++col;
      construct_cdad(i,0,col,1.0);++col;
      construct_adau(i,0,col,1.0);++col;

      //construct singles to form outgoing doubles
      row = 1;

      while(HamOp::ostates[i-1][row].size() == 1){//create up

         construct_cu(i,row,col,1.0);
         ++row;
         ++col;

      }

      row = 1;

      while(HamOp::ostates[i-1][row].size() == 1){//create down with sign

         construct_cds(i,row,col,1.0);
         ++row;
         ++col;

      }

      row = 1;

      while(HamOp::ostates[i-1][row].size() == 1){//annihilate up

         construct_au(i,row,col,1.0);
         ++row;
         ++col;

      }

      row = 1;

      while(HamOp::ostates[i-1][row].size() == 1){//annihilate down with sign

         construct_ads(i,row,col,1.0);
         ++row;
         ++col;

      }

      //copy the doubles from previous tensor: identity
      while(HamOp::ostates[i-1][row].size() == 2){

         construct_id(i,row,col,1.0);
         ++row;
         ++col;

      }

      //HERE STARTS THE COMPLEMENTARY OPERATOR STUFF!
      while(col < HamOp::ostates[i].size() - 1){

         int j = HamOp::ostates[i][col].gsite(0);
         int sj = HamOp::ostates[i][col].gspin(0);
         int aj = HamOp::ostates[i][col].gact(0);

         //first row
         if(sj == 0 && aj == 0)
            construct_tcuf(i,0,col,t(i,j),V(i,i,j,i));
         else if(sj == 1 && aj == 0)
            construct_tcdf(i,0,col,t(i,j),V(i,i,j,i));
         else if(sj == 0 && aj == 1)
            construct_tauf(i,0,col,t(i,j),V(i,i,j,i));
         else
            construct_tadf(i,0,col,t(i,j),V(i,i,j,i));

         //singles coming in:
         row = 1;

         while(HamOp::ostates[i-1][row].size() == 1){

            std::vector<double> val(2);

            std::vector<int> v = Ostate::get_double_complement(i,HamOp::ostates[i-1][row],HamOp::ostates[i][col],V,val);

            if(v.size() == 1){

               if(v[0] == 0)
                  construct_adau(i,row,col,-val[0]);//extra minus sign!
               else if(v[0] == 1)
                  construct_cucd(i,row,col,val[0]);
               else if(v[0] == 2)
                  construct_cuad(i,row,col,val[0]);
               else
                  construct_cdau(i,row,col,val[0]);

            }
            else if(v.size() == 2)//with sign because in the middle!
               construct_pair_s(i,row,col,val);

            ++row;

         }

         //pairs coming in
         while(HamOp::ostates[i-1][row].size() == 2){

            double val;

            std::vector<int> v = Ostate::get_single_complement(i,HamOp::ostates[i-1][row],HamOp::ostates[i][col],V,val);

            if(v.size() > 0){

               if(v[0] == 0 && v[1] == 0)//create spin up
                  construct_cus(i,row,col,val);
               else if(v[0] == 1 && v[1] == 0)//create spin down
                  construct_cd(i,row,col,val);
               else if(v[0] == 0 && v[1] == 1)//annihilate spin up
                  construct_aus(i,row,col,val);
               else//annihilate spin down
                  construct_ad(i,row,col,val);

            }

            ++row;

         }

         //signs: find row and col which are connected
         while(HamOp::ostates[i-1][row] != HamOp::ostates[i][col])
            ++row;

         construct_s(i,row,col,1.0);

         col++;

      }

      //last col! closing down everything!

      //first row
      construct_local(i,0,col,t(i,i),V(i,i,i,i));

      //close down the singles coming in with a triplet
      row = 1;

      while(HamOp::ostates[i-1][row].size() == 1){

         //incoming operator
         int k = HamOp::ostates[i-1][row].gsite(0);
         int sk = HamOp::ostates[i-1][row].gspin(0);
         int ak = HamOp::ostates[i-1][row].gact(0);

         if(sk == 0 && ak == 0)//create up coming in
            construct_tcul(i,row,col,V(k,i,i,i));
         else if(sk == 1 && ak == 0)//create down coming in
            construct_tcdl(i,row,col,-V(i,k,i,i));
         else if(sk == 0 && ak == 1)//annihilate up coming in
            construct_taul(i,row,col,V(i,i,k,i));
         else //annihilate down coming in
            construct_tadl(i,row,col,-V(i,i,i,k));

         ++row;

      }

      //close down the doubles coming in with a pair
      while(HamOp::ostates[i-1][row].size() == 2){

         std::vector<double> val(2);

         std::vector<int> v = Ostate::get_closing_pair_in(i,HamOp::ostates[i-1][row],V,val);

         if(v.size() == 1){

            if(v[0] == 0)
               construct_adau(i,row,col,-val[0]);//extra minus sign!
            else if(v[0] == 1)
               construct_cucd(i,row,col,val[0]);
            else if(v[0] == 2)
               construct_cuad(i,row,col,val[0]);
            else
               construct_cdau(i,row,col,val[0]);

         }
         else if(v.size() == 2)
            construct_pair(i,row,col,val);

         ++row;

      }

      //close down the complementary operators of this site
      //basically the first 4 incoming states should be closed down
      construct_cu(i,row,col,1.0);
      ++row;

      construct_cds(i,row,col,1.0);
      ++row;

      construct_au(i,row,col,1.0);
      ++row;

      construct_ads(i,row,col,1.0);
      ++row;

      //finally construct the identity on the lower right element
      construct_id(i,HamOp::ostates[i-1].size() - 1,HamOp::ostates[i].size() - 1,1.0);

   }

   //first row completely closed
   construct_local(L/2,0,0,t(L/2,L/2),V(L/2,L/2,L/2,L/2));

   //construct complementary triples:
   row = 1;

   while(HamOp::ostates[L/2 - 1][row].size() == 1){

      //incoming operator
      int k = HamOp::ostates[L/2 - 1][row].gsite(0);
      int sk = HamOp::ostates[L/2 - 1][row].gspin(0);
      int ak = HamOp::ostates[L/2 - 1][row].gact(0);

      if(sk == 0 && ak == 0)//create up coming in
         construct_tcul(L/2,row,0,V(k,L/2,L/2,L/2));
      else if(sk == 1 && ak == 0)//create down coming in
         construct_tcdl(L/2,row,0,-V(L/2,k,L/2,L/2));
      else if(sk == 0 && ak == 1)//annihilate up coming in
         construct_taul(L/2,row,0,V(L/2,L/2,k,L/2));
      else //annihilate down coming in
         construct_tadl(L/2,row,0,-V(L/2,L/2,L/2,k));

      ++row;

   }

   //close down the doubles coming in with a pair
   while(HamOp::ostates[L/2 - 1][row].size() == 2){

      std::vector<double> val(2);

      std::vector<int> v = Ostate::get_closing_pair_in(L/2,HamOp::ostates[L/2 - 1][row],V,val);

      if(v.size() == 1){

         if(v[0] == 0)
            construct_adau(L/2,row,0,-val[0]);//extra minus sign!
         else if(v[0] == 1)
            construct_cucd(L/2,row,0,val[0]);
         else if(v[0] == 2)
            construct_cuad(L/2,row,0,val[0]);
         else
            construct_cdau(L/2,row,0,val[0]);

      }
      else if(v.size() == 2)
         construct_pair(L/2,row,0,val);

      ++row;

   }

   //close down the complementary triple operators of this site
   //basically the first 4 incoming states should be closed down
   construct_cu(L/2,row,0,1.0);
   ++row;

   construct_cds(L/2,row,0,1.0);
   ++row;

   construct_au(L/2,row,0,1.0);
   ++row;

   construct_ads(L/2,row,0,1.0);
   ++row;

   //finally construct the identity for the already closed terms coming in
   construct_id(L/2,HamOp::ostates[L/2 - 1].size() - 1,0,1.0);

   //next cols: outgoing singles
   col = 1;

   while(HamOp::ostates[L/2][col].size() == 1){

      int j = HamOp::ostates[L/2][col].gsite(0);
      int sj = HamOp::ostates[L/2][col].gspin(0);
      int aj = HamOp::ostates[L/2][col].gact(0);

      //first row
      if(sj == 0 && aj == 0)
         construct_tcuf(L/2,0,col,t(L/2,j),V(L/2,L/2,j,L/2));
      else if(sj == 1 && aj == 0)
         construct_tcdf(L/2,0,col,t(L/2,j),V(L/2,L/2,j,L/2));
      else if(sj == 0 && aj == 1)
         construct_tauf(L/2,0,col,t(L/2,j),V(L/2,L/2,j,L/2));
      else
         construct_tadf(L/2,0,col,t(L/2,j),V(L/2,L/2,j,L/2));

      //singles coming in:
      row = 1;

      while(HamOp::ostates[L/2 - 1][row].size() == 1){

         std::vector<double> val(2);

         std::vector<int> v = Ostate::get_double_complement(L/2,HamOp::ostates[L/2 - 1][row],HamOp::ostates[L/2][col],V,val);

         if(v.size() == 1){

            if(v[0] == 0)
               construct_adau(L/2,row,col,-val[0]);//extra minus sign!
            else if(v[0] == 1)
               construct_cucd(L/2,row,col,val[0]);
            else if(v[0] == 2)
               construct_cuad(L/2,row,col,val[0]);
            else
               construct_cdau(L/2,row,col,val[0]);

         }
         else if(v.size() == 2)//with sign because in the middle!
            construct_pair_s(L/2,row,col,val);

         ++row;

      }

      //pairs coming in
      while(HamOp::ostates[L/2 - 1][row].size() == 2){

         double val;

         std::vector<int> v = Ostate::get_single_complement(L/2,HamOp::ostates[L/2 - 1][row],HamOp::ostates[L/2][col],V,val);

         if(v.size() > 0){

            if(v[0] == 0 && v[1] == 0)//create spin up
               construct_cus(L/2,row,col,val);
            else if(v[0] == 1 && v[1] == 0)//create spin down
               construct_cd(L/2,row,col,val);
            else if(v[0] == 0 && v[1] == 1)//annihilate spin up
               construct_aus(L/2,row,col,val);
            else//annihilate spin down
               construct_ad(L/2,row,col,val);

         }

         ++row;

      }

      //signs: find row and col which are connected
      while(HamOp::ostates[L/2 - 1][row] != HamOp::ostates[L/2][col])
         ++row;

      construct_s(L/2,row,col,1.0);

      col++;

   }

   //now outgoing pairs!
   while(HamOp::ostates[L/2][col].size() == 2){

      int k = HamOp::ostates[L/2][col].gsite(0);
      int sk = HamOp::ostates[L/2][col].gspin(0);
      int ak = HamOp::ostates[L/2][col].gact(0);

      int l = HamOp::ostates[L/2][col].gsite(1);
      int sl = HamOp::ostates[L/2][col].gspin(1);
      int al = HamOp::ostates[L/2][col].gact(1);

      //first row: id coming in, construct pairs
      std::vector<double> val(2);

      std::vector<int> v = Ostate::get_closing_pair_out(L/2,HamOp::ostates[L/2][col],V,val);

      if(v.size() == 1){

         if(v[0] == 0)
            construct_adau(L/2,0,col,-val[0]);//extra minus sign!
         else if(v[0] == 1)
            construct_cucd(L/2,0,col,val[0]);
         else if(v[0] == 2)
            construct_cuad(L/2,0,col,val[0]);
         else
            construct_cdau(L/2,0,col,val[0]);

      }
      else if(v.size() == 2)
         construct_pair(L/2,0,col,val);

      //singles coming in: construct single to match outgoing pair
      row = 1;

      while(HamOp::ostates[L/2 - 1][row].size() == 1){

         double val;

         std::vector<int> v = Ostate::get_single_complement(L/2,HamOp::ostates[L/2][col],HamOp::ostates[L/2 - 1][row],V,val);

         if(v.size() > 0){

            if(v[0] == 0 && v[1] == 0)//create spin up
               construct_cu(L/2,row,col,val);
            else if(v[0] == 1 && v[1] == 0)//create spin down
               construct_cds(L/2,row,col,val);
            else if(v[0] == 0 && v[1] == 1)//annihilate spin up
               construct_au(L/2,row,col,val);
            else//annihilate spin down
               construct_ads(L/2,row,col,val);

         }

         ++row;

      }

      //pairs coming in: convert incoming to outgoing pairs

      while(HamOp::ostates[L/2 - 1][row].size() == 2){

         double val;

         int op = Ostate::transfer_pair_in_pair_out(HamOp::ostates[L/2 - 1][row],HamOp::ostates[L/2][col],V,val);

         if(op == 1)
            construct_id(L/2,row,col,val);

         ++row;

      }

      ++col;

   }

   //then the outgoing "incoming" singles
   while(col < HamOp::ostates[L/2].size() - 1){

      int ci = HamOp::ostates[L/2][col].gsite(0);
      int cs = HamOp::ostates[L/2][col].gspin(0);
      int ca = HamOp::ostates[L/2][col].gact(0);

      //construct single operators on the first row
      if(ci == L/2){

         if(cs == 0 && ca == 0)
            construct_cus(L/2,0,col,1.0);
         else if(cs == 1 && ca == 0)
            construct_cd(L/2,0,col,1.0);
         else if(cs == 0 && ca == 1)
            construct_aus(L/2,0,col,1.0);
         else
            construct_ad(L/2,0,col,1.0);

      }

      row = 1;

      while(HamOp::ostates[L/2 - 1][row].size() == 1){

         int ri = HamOp::ostates[L/2 - 1][row].gsite(0);
         int rs = HamOp::ostates[L/2 - 1][row].gspin(0);
         int ra = HamOp::ostates[L/2 - 1][row].gact(0);

         if(ri == ci && rs == cs && ra == ca)
            construct_s(L/2,row,col,1.0);

         ++row;

      }

      ++col;

   }

   //finally construct id for incoming id
   construct_id(L/2,0,HamOp::ostates[L/2].size() - 1,1.0);

   //the rest of the blocks: L/2 + 1 --> L - 1
   for(int i = L/2 + 1;i < L - 1;++i){

      col = 0;

      //first col is closed
      construct_id(i,0,0,1.0);

      //for incoming complementary triples, close down with single
      row = 1;

      while(HamOp::ostates[i - 1][row].size() == 1){

         int ri = HamOp::ostates[i-1][row].gsite(0);
         int rs = HamOp::ostates[i-1][row].gspin(0);
         int ra = HamOp::ostates[i-1][row].gact(0);

         if(ri == i){

            if(rs == 0 && ra == 0)
               construct_cu(i,row,col,1.0);
            else if(rs == 1 && ra == 0)
               construct_cds(i,row,col,1.0);
            else if(rs == 0 && ra == 1)
               construct_au(i,row,col,1.0);
            else
               construct_ads(i,row,col,1.0);

         }

         ++row;

      }

      //close down incoming complementary doubles with pairs
      while(HamOp::ostates[i - 1][row].size() == 2){

         int ri_1 = HamOp::ostates[i-1][row].gsite(0);
         int rs_1 = HamOp::ostates[i-1][row].gspin(0);
         int ra_1 = HamOp::ostates[i-1][row].gact(0);

         int ri_2 = HamOp::ostates[i-1][row].gsite(1);
         int rs_2 = HamOp::ostates[i-1][row].gspin(1);
         int ra_2 = HamOp::ostates[i-1][row].gact(1);

         if(ri_1 == i && ri_2 == i){

            if(ra_1 == 0){//first crea

               if(rs_1 == 0){//first create up

                  if(rs_2 == 1 && ra_2 == 0)
                     construct_cucd(i,row,col,1.0);
                  else if(rs_2 == 0 && ra_2 == 1)
                     construct_cuau(i,row,col,1.0);
                  else if(rs_2 == 1 && ra_2 == 1)
                     construct_cuad(i,row,col,1.0);

               }
               else{//first create down

                  if(rs_2 == 0 && ra_2 == 1)
                     construct_cdau(i,row,col,1.0);
                  else if(rs_2 == 1 && ra_2 == 1)
                     construct_cdad(i,row,col,1.0);

               }

            }
            else{//first anni

               if(rs_1 == 1){//first anni down

                  if(rs_2 == 0 && ra_2 == 1)//second anni up: extra minus sign!
                     construct_adau(i,row,col,1.0);

               }

            }

         }

         ++row;

      }

      //incoming singles, close with triples!
      while(row < HamOp::ostates[i - 1].size() - 1){

         //incoming operator
         int k = HamOp::ostates[i-1][row].gsite(0);
         int sk = HamOp::ostates[i-1][row].gspin(0);
         int ak = HamOp::ostates[i-1][row].gact(0);

         if(sk == 0 && ak == 0)//create up coming in
            construct_tcul(i,row,0,V(k,i,i,i));
         else if(sk == 1 && ak == 0)//create down coming in
            construct_tcdl(i,row,0,-V(i,k,i,i));
         else if(sk == 0 && ak == 1)//annihilate up coming in
            construct_taul(i,row,0,V(i,i,k,i));
         else //annihilate down coming in
            construct_tadl(i,row,0,-V(i,i,i,k));

         ++row;

      }

      //construct local terms
      construct_local(i,HamOp::ostates[i - 1].size() - 1,0,t(i,i),V(i,i,i,i));

      //outgoing complementary triples
      col = 1;

      while(HamOp::ostates[i][col].size() == 1){

         int ci = HamOp::ostates[i][col].gsite(0);
         int cs = HamOp::ostates[i][col].gspin(0);
         int ca = HamOp::ostates[i][col].gact(0);

         //incoming complementary triples
         row = 1;

         while(HamOp::ostates[i-1][row].size() == 1){

            int ri = HamOp::ostates[i-1][row].gsite(0);
            int rs = HamOp::ostates[i-1][row].gspin(0);
            int ra = HamOp::ostates[i-1][row].gact(0);

            if(ci == ri && cs == rs && ca == ra)
               construct_s(i,row,col,1.0);

            ++row;

         }

         //incoming complementary pairs: add single
         while(HamOp::ostates[i-1][row].size() == 2){

            int ri1 = HamOp::ostates[i-1][row].gsite(0);
            int rs1 = HamOp::ostates[i-1][row].gspin(0);
            int ra1 = HamOp::ostates[i-1][row].gact(0);

            int ri2 = HamOp::ostates[i-1][row].gsite(1);
            int rs2 = HamOp::ostates[i-1][row].gspin(1);
            int ra2 = HamOp::ostates[i-1][row].gact(1);

            if(ri1 == i){//if first site of the complementary is this site i:

               if(ri2 == ci && rs2 == cs && ra2 == ca){//and the second matches up with the outgoing single

                  if(rs1 == 0 && ra1 == 0)
                     construct_cus(i,row,col,1.0);
                  else if(rs1 == 1 && ra1 == 0)
                     construct_cd(i,row,col,1.0);
                  else if(rs1 == 0 && ra1 == 1)
                     construct_aus(i,row,col,1.0);
                  else
                     construct_ad(i,row,col,1.0);

               }

            }

            ++row;

         }

         //singles coming in:
         while(row < HamOp::ostates[i-1].size() - 1){

            std::vector<double> val(2);

            std::vector<int> v = Ostate::get_double_complement(i,HamOp::ostates[i-1][row],HamOp::ostates[i][col],V,val);

            if(v.size() == 1){

               if(v[0] == 0)
                  construct_adau(i,row,col,-val[0]);//extra minus sign!
               else if(v[0] == 1)
                  construct_cucd(i,row,col,val[0]);
               else if(v[0] == 2)
                  construct_cuad(i,row,col,val[0]);
               else
                  construct_cdau(i,row,col,val[0]);

            }
            else if(v.size() == 2)//with sign because in the middle!
               construct_pair_s(i,row,col,val);

            ++row;

         }

         if(cs == 0 && ca == 0)
            construct_tcuf(i,row,col,t(i,ci),V(i,i,ci,i));
         else if(cs == 1 && ca == 0)
            construct_tcdf(i,row,col,t(i,ci),V(i,i,ci,i));
         else if(cs == 0 && ca == 1)
            construct_tauf(i,row,col,t(i,ci),V(i,i,ci,i));
         else
            construct_tadf(i,row,col,t(i,ci),V(i,i,ci,i));

         ++col;

      }

      while(HamOp::ostates[i][col].size() == 2){

         //now start from the bottom, with the last row
         row = HamOp::ostates[i-1].size() - 1;

         //last row: id coming in, construct pairs
         std::vector<double> val(2);

         std::vector<int> v = Ostate::get_closing_pair_out(i,HamOp::ostates[i][col],V,val);

         if(v.size() == 1){

            if(v[0] == 0)
               construct_adau(i,row,col,-val[0]);//extra minus sign!
            else if(v[0] == 1)
               construct_cucd(i,row,col,val[0]);
            else if(v[0] == 2)
               construct_cuad(i,row,col,val[0]);
            else
               construct_cdau(i,row,col,val[0]);

         }
         else if(v.size() == 2)
            construct_pair(i,row,col,val);

         row--;

         //now decrease the row: incoming singles
         while(HamOp::ostates[i-1][row].size() == 1){

            double val;

            std::vector<int> v = Ostate::get_single_complement(i,HamOp::ostates[i][col],HamOp::ostates[i-1][row],V,val);

            if(v.size() > 0){

               if(v[0] == 0 && v[1] == 0)//create spin up
                  construct_cu(i,row,col,val);
               else if(v[0] == 1 && v[1] == 0)//create spin down
                  construct_cds(i,row,col,val);
               else if(v[0] == 0 && v[1] == 1)//annihilate spin up
                  construct_au(i,row,col,val);
               else//annihilate spin down
                  construct_ads(i,row,col,val);

            }

            --row;

         }

         //incoming complementary pairs: put in id on the right place
         while(HamOp::ostates[i-1][row].size() == 2){

            if(HamOp::ostates[i-1][row] == HamOp::ostates[i][col])
               construct_id(i,row,col,1.0);

            --row;

         }

         ++col;

      }

      //outgoing singles
      while(col < HamOp::ostates[i].size() - 1){

         int ci = HamOp::ostates[i][col].gsite(0);
         int cs = HamOp::ostates[i][col].gspin(0);
         int ca = HamOp::ostates[i][col].gact(0);

         row = HamOp::ostates[i - 1].size() - 1;

         //construct correct operator:
         if(ci == i){

            if(cs == 0 && ca == 0)
               construct_cus(i,row,col,1.0);
            else if(cs == 1 && ca == 0)
               construct_cd(i,row,col,1.0);
            else if(cs == 0 && ca == 1)
               construct_aus(i,row,col,1.0);
            else
               construct_ad(i,row,col,1.0);

         }

         row--;

         //construct sign for incoming singles to transfer to outgoing singles
         while(HamOp::ostates[i - 1][row].size() == 1){

            if(HamOp::ostates[i - 1][row] == HamOp::ostates[i][col])
               construct_s(i,row,col,1.0);

            --row;

         }

         ++col;

      }

      //finally construct id in bottom right
      construct_id(i,HamOp::ostates[i - 1].size() - 1,HamOp::ostates[i].size() - 1,1.0);

   }

   //finally the last mpo:
   construct_id(L-1,0,0,1.0);

   //for incoming complementary triples, close down with single
   row = 1;

   while(HamOp::ostates[L - 2][row].size() == 1){

      int ri = HamOp::ostates[L-2][row].gsite(0);
      int rs = HamOp::ostates[L-2][row].gspin(0);
      int ra = HamOp::ostates[L-2][row].gact(0);

      if(ri == L - 1){

         if(rs == 0 && ra == 0)
            construct_cu(L-1,row,0,1.0);
         else if(rs == 1 && ra == 0)
            construct_cds(L-1,row,0,1.0);
         else if(rs == 0 && ra == 1)
            construct_au(L-1,row,0,1.0);
         else
            construct_ads(L-1,row,0,1.0);

      }

      ++row;

   }

   //close down incoming complementary doubles with pairs
   while(HamOp::ostates[L-2][row].size() == 2){

      int ri_1 = HamOp::ostates[L-2][row].gsite(0);
      int rs_1 = HamOp::ostates[L-2][row].gspin(0);
      int ra_1 = HamOp::ostates[L-2][row].gact(0);

      int ri_2 = HamOp::ostates[L-2][row].gsite(1);
      int rs_2 = HamOp::ostates[L-2][row].gspin(1);
      int ra_2 = HamOp::ostates[L-2][row].gact(1);

      if(ri_1 == L - 1 && ri_2 == L - 1){

         if(ra_1 == 0){//first crea

            if(rs_1 == 0){//first create up

               if(rs_2 == 1 && ra_2 == 0)
                  construct_cucd(L-1,row,0,1.0);
               else if(rs_2 == 0 && ra_2 == 1)
                  construct_cuau(L-1,row,0,1.0);
               else if(rs_2 == 1 && ra_2 == 1)
                  construct_cuad(L-1,row,0,1.0);

            }
            else{//first create down

               if(rs_2 == 0 && ra_2 == 1)
                  construct_cdau(L-1,row,0,1.0);
               else if(rs_2 == 1 && ra_2 == 1)
                  construct_cdad(L-1,row,0,1.0);

            }

         }
         else{//first anni

            if(rs_1 == 1){//first anni down

               if(rs_2 == 0 && ra_2 == 1)//second anni up: extra minus sign!
                  construct_adau(L-1,row,0,1.0);

            }

         }

      }

      ++row;

   }

   //incoming singles, close with triples!
   while(row < HamOp::ostates[L-2].size() - 1){

      //incoming operator
      int k = HamOp::ostates[L-2][row].gsite(0);
      int sk = HamOp::ostates[L-2][row].gspin(0);
      int ak = HamOp::ostates[L-2][row].gact(0);

      if(sk == 0 && ak == 0)//create up coming in
         construct_tcul(L-1,row,0,V(k,L - 1,L - 1,L - 1));
      else if(sk == 1 && ak == 0)//create down coming in
         construct_tcdl(L-1,row,0,-V(L - 1,k,L - 1,L - 1));
      else if(sk == 0 && ak == 1)//annihilate up coming in
         construct_taul(L-1,row,0,V(L - 1,L - 1,k,L - 1));
      else //annihilate down coming in
         construct_tadl(L-1,row,0,-V(L - 1,L - 1,L - 1,k));

      ++row;

   }

   //construct local terms
   construct_local(L-1,HamOp::ostates[L-2].size() - 1,0,t(L - 1,L - 1),V(L - 1,L - 1,L - 1,L - 1));

}

/**
 * delete all the lists and stuff
 */
void Operator::clear(){

   for(int i = 0;i < L;++i)
      for(int row = 0;row < dim[i][0];++row)
         for(int col = 0;col < dim[i][1];++col)
            if(sparse[i][row*dim[i][1] + col])
               delete op[i][row*dim[i][1] + col];

   for(int i = 0;i < L;++i)
      delete [] op[i];

   for(int i = 0;i < L;++i)
      delete [] sparse[i];

}

/**
 * construct identity operator in mpo O
 */
void Operator::construct_id(int site,int row,int col,double val){

   op[site][row*dim[site][1] + col] = new QSDArray<2>(Quantum::zero(),make_array(qp,-qp));
   sparse[site][row*dim[site][1] + col] = true;

   DArray<2> Ip(1,1);
   Ip = val;

   op[site][row*dim[site][1] + col]->insert(shape(0,0),Ip);
   op[site][row*dim[site][1] + col]->insert(shape(1,1),Ip);
   op[site][row*dim[site][1] + col]->insert(shape(2,2),Ip);
   op[site][row*dim[site][1] + col]->insert(shape(3,3),Ip);

}
/**
 * construct identity operator with fermion sign
 */
void Operator::construct_s(int site,int row,int col,double val){

   op[site][row*dim[site][1] + col] = new QSDArray<2>(Quantum::zero(),make_array(qp,-qp));
   sparse[site][row*dim[site][1] + col] = true;

   DArray<2> Ip(1,1);

   Ip = val;
   op[site][row*dim[site][1] + col]->insert(shape(0,0),Ip);
   op[site][row*dim[site][1] + col]->insert(shape(3,3),Ip);

   Ip = -val;
   op[site][row*dim[site][1] + col]->insert(shape(1,1),Ip);
   op[site][row*dim[site][1] + col]->insert(shape(2,2),Ip);

}

/**
 * construct creator of up spin
 */
void Operator::construct_cu(int site,int row,int col,double val){

   op[site][row*dim[site][1] + col] = new QSDArray<2>(Quantum(1,0),make_array(qp,-qp));
   sparse[site][row*dim[site][1] + col] = true;

   DArray<2> Ip(1,1);
   Ip = val;

   op[site][row*dim[site][1] + col]->insert(shape(1,0),Ip);
   op[site][row*dim[site][1] + col]->insert(shape(3,2),Ip);

}

/**
 * construct creator of up spin with sign down
 */
void Operator::construct_cus(int site,int row,int col,double val){

   op[site][row*dim[site][1] + col] = new QSDArray<2>(Quantum(1,0),make_array(qp,-qp));
   sparse[site][row*dim[site][1] + col] = true;

   DArray<2> Ip(1,1);
   Ip = val;

   //a^+_up (-1)^n_down
   op[site][row*dim[site][1] + col]->insert(shape(1,0),Ip);

   Ip = -val;

   op[site][row*dim[site][1] + col]->insert(shape(3,2),Ip);

}

/**
 * construct creator of down spin
 */
void Operator::construct_cd(int site,int row,int col,double val){

   op[site][row*dim[site][1] + col] = new QSDArray<2>(Quantum(0,1),make_array(qp,-qp));
   sparse[site][row*dim[site][1] + col] = true;

   DArray<2> Ip(1,1);
   Ip = val;

   //a^dagger_down 
   op[site][row*dim[site][1] + col]->insert(shape(2,0),Ip);
   op[site][row*dim[site][1] + col]->insert(shape(3,1),Ip);

}

/**
 * construct creator of down spin with up spin sign
 */
void Operator::construct_cds(int site,int row,int col,double val){

   op[site][row*dim[site][1] + col] = new QSDArray<2>(Quantum(0,1),make_array(qp,-qp));
   sparse[site][row*dim[site][1] + col] = true;

   DArray<2> Ip(1,1);

   //a^dagger_down 
   Ip = val;
   op[site][row*dim[site][1] + col]->insert(shape(2,0),Ip);

   Ip = -val;
   op[site][row*dim[site][1] + col]->insert(shape(3,1),Ip);

}

/**
 * construct annihilator of up spin
 */
void Operator::construct_au(int site,int row,int col,double val){

   op[site][row*dim[site][1] + col] = new QSDArray<2>(Quantum(-1,0),make_array(qp,-qp));
   sparse[site][row*dim[site][1] + col] = true;

   DArray<2> Ip(1,1);
   Ip = val;

   //a_up (-1)^n_down
   op[site][row*dim[site][1] + col]->insert(shape(0,1),Ip);
   op[site][row*dim[site][1] + col]->insert(shape(2,3),Ip);

}

/**
 * construct annihilator of up spin with down sign
 */
void Operator::construct_aus(int site,int row,int col,double val){

   op[site][row*dim[site][1] + col] = new QSDArray<2>(Quantum(-1,0),make_array(qp,-qp));
   sparse[site][row*dim[site][1] + col] = true;

   DArray<2> Ip(1,1);

   //a_up (-1)^n_down
   Ip = val;
   op[site][row*dim[site][1] + col]->insert(shape(0,1),Ip);

   Ip = -val;
   op[site][row*dim[site][1] + col]->insert(shape(2,3),Ip);

}

/**
 * construct annihilator of down spin
 */
void Operator::construct_ad(int site,int row,int col,double val){

   op[site][row*dim[site][1] + col] = new QSDArray<2>(Quantum(0,-1),make_array(qp,-qp));
   sparse[site][row*dim[site][1] + col] = true;

   DArray<2> Ip(1,1);
   Ip = val;

   //a_down
   op[site][row*dim[site][1] + col]->insert(shape(0,2),Ip);
   op[site][row*dim[site][1] + col]->insert(shape(1,3),Ip);

}

/**
 * construct annihilator of down spin with sign for up spin
 */
void Operator::construct_ads(int site,int row,int col,double val){

   op[site][row*dim[site][1] + col] = new QSDArray<2>(Quantum(0,-1),make_array(qp,-qp));
   sparse[site][row*dim[site][1] + col] = true;

   DArray<2> Ip(1,1);
   Ip = val;

   //a_down
   Ip = val;
   op[site][row*dim[site][1] + col]->insert(shape(0,2),Ip);

   Ip = -val;
   op[site][row*dim[site][1] + col]->insert(shape(1,3),Ip);

}

/**
 * construct creator of an up down pair on site
 */
void Operator::construct_cucd(int site,int row,int col,double val){

   op[site][row*dim[site][1] + col] = new QSDArray<2>(Quantum(1,1),make_array(qp,-qp));
   sparse[site][row*dim[site][1] + col] = true;

   DArray<2> Ip(1,1);
   Ip = val;

   //a_down
   op[site][row*dim[site][1] + col]->insert(shape(3,0),Ip);

}

/**
 * construct n_up operator
 */
void Operator::construct_cuau(int site,int row,int col,double val){

   op[site][row*dim[site][1] + col] = new QSDArray<2>(Quantum::zero(),make_array(qp,-qp));
   sparse[site][row*dim[site][1] + col] = true;

   DArray<2> Ip(1,1);
   Ip = val;

   op[site][row*dim[site][1] + col]->insert(shape(1,1),Ip);
   op[site][row*dim[site][1] + col]->insert(shape(3,3),Ip);

}

/**
 * construct create up annihilate down
 */
 
void Operator::construct_cuad(int site,int row,int col,double val){

   op[site][row*dim[site][1] + col] = new QSDArray<2>(Quantum(1,-1),make_array(qp,-qp));
   sparse[site][row*dim[site][1] + col] = true;

   DArray<2> Ip(1,1);
   Ip = val;

   op[site][row*dim[site][1] + col]->insert(shape(1,2),Ip);

}

/**
 * construct create down annihilate up
 */
void Operator::construct_cdau(int site,int row,int col,double val){

   op[site][row*dim[site][1] + col] = new QSDArray<2>(Quantum(-1,1),make_array(qp,-qp));
   sparse[site][row*dim[site][1] + col] = true;

   DArray<2> Ip(1,1);
   Ip = val;

   op[site][row*dim[site][1] + col]->insert(shape(2,1),Ip);

}

/**
 * construct create down annihilate down --> n_down
 */
void Operator::construct_cdad(int site,int row,int col,double val){

   op[site][row*dim[site][1] + col] = new QSDArray<2>(Quantum::zero(),make_array(qp,-qp));
   sparse[site][row*dim[site][1] + col] = true;

   DArray<2> Ip(1,1);
   Ip = val;

   op[site][row*dim[site][1] + col]->insert(shape(2,2),Ip);
   op[site][row*dim[site][1] + col]->insert(shape(3,3),Ip);

}

/**
 * construct annihilate up annihilate down: extra minus sign!
 */
void Operator::construct_adau(int site,int row,int col,double val){

   op[site][row*dim[site][1] + col] = new QSDArray<2>(Quantum(-1,-1),make_array(qp,-qp));
   sparse[site][row*dim[site][1] + col] = true;

   DArray<2> Ip(1,1);
   Ip = -val;

   op[site][row*dim[site][1] + col]->insert(shape(0,3),Ip);

}

/**
 * construct complementary operator for a_up: behaves as a creator of up
 */
void Operator::construct_tauf(int site,int row,int col,double tval,double Vval){

   op[site][row*dim[site][1] + col] = new QSDArray<2>(Quantum(1,0),make_array(qp,-qp));
   sparse[site][row*dim[site][1] + col] = true;

   DArray<2> Ip(1,1);

   Ip = tval;
   op[site][row*dim[site][1] + col]->insert(shape(1,0),Ip);

   Ip = -tval - Vval;
   op[site][row*dim[site][1] + col]->insert(shape(3,2),Ip);

}

/**
 * construct complementary operator for a_up: behaves as a creator of up
 */
void Operator::construct_taul(int site,int row,int col,double Vval){

   op[site][row*dim[site][1] + col] = new QSDArray<2>(Quantum(1,0),make_array(qp,-qp));
   sparse[site][row*dim[site][1] + col] = true;

   DArray<2> Ip(1,1);

   Ip = Vval;
   op[site][row*dim[site][1] + col]->insert(shape(3,2),Ip);

}

/**
 * construct complementary operator for a_down: behaves as a creator of down 
 */
void Operator::construct_tadf(int site,int row,int col,double tval,double Vval){

   op[site][row*dim[site][1] + col] = new QSDArray<2>(Quantum(0,1),make_array(qp,-qp));
   sparse[site][row*dim[site][1] + col] = true;

   DArray<2> Ip(1,1);

   Ip = tval;
   op[site][row*dim[site][1] + col]->insert(shape(2,0),Ip);

   Ip = tval + Vval;
   op[site][row*dim[site][1] + col]->insert(shape(3,1),Ip);

}

/**
 * construct complementary operator for a_down: behaves as a creator of down 
 */
void Operator::construct_tadl(int site,int row,int col,double Vval){

   op[site][row*dim[site][1] + col] = new QSDArray<2>(Quantum(0,1),make_array(qp,-qp));
   sparse[site][row*dim[site][1] + col] = true;

   DArray<2> Ip(1,1);

   Ip = Vval;
   op[site][row*dim[site][1] + col]->insert(shape(3,1),Ip);

}

/**
 * construct complementary operator for a^+_down: behaves as an annihilator of down 
 */
void Operator::construct_tcdf(int site,int row,int col,double tval,double Vval){

   op[site][row*dim[site][1] + col] = new QSDArray<2>(Quantum(0,-1),make_array(qp,-qp));
   sparse[site][row*dim[site][1] + col] = true;

   DArray<2> Ip(1,1);

   Ip = tval;
   op[site][row*dim[site][1] + col]->insert(shape(0,2),Ip);

   Ip = tval + Vval;
   op[site][row*dim[site][1] + col]->insert(shape(1,3),Ip);

}

/**
 * construct complementary operator for a^+_down: behaves as an annihilator of down 
 */
void Operator::construct_tcdl(int site,int row,int col,double Vval){

   op[site][row*dim[site][1] + col] = new QSDArray<2>(Quantum(0,-1),make_array(qp,-qp));
   sparse[site][row*dim[site][1] + col] = true;

   DArray<2> Ip(1,1);

   Ip = Vval;
   op[site][row*dim[site][1] + col]->insert(shape(1,3),Ip);

}

/**
 * construct complementary operator for a^+_up: behaves as an annihilator of up
 */
void Operator::construct_tcuf(int site,int row,int col,double tval,double Vval){

   op[site][row*dim[site][1] + col] = new QSDArray<2>(Quantum(-1,0),make_array(qp,-qp));
   sparse[site][row*dim[site][1] + col] = true;

   DArray<2> Ip(1,1);

   Ip = tval;
   op[site][row*dim[site][1] + col]->insert(shape(0,1),Ip);

   Ip = -tval - Vval;
   op[site][row*dim[site][1] + col]->insert(shape(2,3),Ip);

}

/**
 * construct complementary operator for a^+_up: behaves as an annihilator of up
 */
void Operator::construct_tcul(int site,int row,int col,double Vval){

   op[site][row*dim[site][1] + col] = new QSDArray<2>(Quantum(-1,0),make_array(qp,-qp));
   sparse[site][row*dim[site][1] + col] = true;

   DArray<2> Ip(1,1);

   Ip = Vval;
   op[site][row*dim[site][1] + col]->insert(shape(2,3),Ip);

}

/**
 * construct local term: t(i,i) + V(i,i,i,i)
 */
void Operator::construct_local(int site,int row,int col,double tval,double Vval){

   op[site][row*dim[site][1] + col] = new QSDArray<2>(Quantum::zero(),make_array(qp,-qp));
   sparse[site][row*dim[site][1] + col] = true;

   DArray<2> Ip(1,1);

   Ip = tval;

   op[site][row*dim[site][1] + col]->insert(shape(1,1),Ip);
   op[site][row*dim[site][1] + col]->insert(shape(2,2),Ip);

   Ip = 2*tval + Vval;

   op[site][row*dim[site][1] + col]->insert(shape(3,3),Ip);

}

/**
 * construct pair val[0] a^+_up a_up (-1)^n_down + val[1] a^+_down a_down (-1)^n_up
 */
void Operator::construct_pair_s(int site,int row,int col,const std::vector<double> &val){

   op[site][row*dim[site][1] + col] = new QSDArray<2>(Quantum::zero(),make_array(qp,-qp));
   sparse[site][row*dim[site][1] + col] = true;

   DArray<2> Ip(1,1);

   //up up
   Ip = val[0];
   op[site][row*dim[site][1] + col]->insert(shape(1,1),Ip);

   //down down
   Ip = val[1];
   op[site][row*dim[site][1] + col]->insert(shape(2,2),Ip);

   //both
   Ip = -val[0] - val[1];
   op[site][row*dim[site][1] + col]->insert(shape(3,3),Ip);

}

/**
 * construct pair val[0] a^+_up a_up (-1)^n_down + val[1] a^+_down a_down (-1)^n_up
 */
void Operator::construct_pair(int site,int row,int col,const std::vector<double> &val){

   op[site][row*dim[site][1] + col] = new QSDArray<2>(Quantum::zero(),make_array(qp,-qp));
   sparse[site][row*dim[site][1] + col] = true;

   DArray<2> Ip(1,1);

   //up up
   Ip = val[0];
   op[site][row*dim[site][1] + col]->insert(shape(1,1),Ip);

   //down down
   Ip = val[1];
   op[site][row*dim[site][1] + col]->insert(shape(2,2),Ip);

   //both
   Ip = val[0] + val[1];
   op[site][row*dim[site][1] + col]->insert(shape(3,3),Ip);

}

/**
 * print all the operators, as a test thingy
 */
void Operator::print(){

   cout << endl;
   cout << "******************************************************************************" << endl;
   cout << endl;
   cout << "Operators on site " << 0 << endl;
   cout << endl;
   cout << "******************************************************************************" << endl;
   cout << endl;

   for(int col = 0;col < dim[0][1];++col){

         if(sparse[0][col]){

            cout << endl;
            cout << 0 << "\t" << col << endl;
            cout << endl;
            cout << HamOp::ostates[0][0] << endl;
            cout << HamOp::ostates[0][col] << endl;
            cout << endl;
            cout << *op[0][col] << endl;
            cout << endl;

         }

   }

   for(int i = 1;i < L - 1;++i){

      cout << endl;
      cout << "******************************************************************************" << endl;
      cout << endl;
      cout << "Operators on site " << i << endl;
      cout << endl;
      cout << "******************************************************************************" << endl;
      cout << endl;

      for(int row = 0;row < dim[i][0];++row)
         for(int col = 0;col < dim[i][1];++col){

            if(sparse[i][row*dim[i][1] + col]){

               cout << endl;
               cout << row << "\t" << col << endl;
               cout << endl;
               cout << HamOp::ostates[i - 1][row] << endl;
               cout << HamOp::ostates[i][col] << endl;
               cout << endl;
               cout << *op[i][row*dim[i][1] + col] << endl;
               cout << endl;

            }

         }

   }

   cout << endl;
   cout << "******************************************************************************" << endl;
   cout << endl;
   cout << "Operators on site " << L - 1 << endl;
   cout << endl;
   cout << "******************************************************************************" << endl;
   cout << endl;

   for(int row = 0;row < dim[L - 1][0];++row)
      for(int col = 0;col < dim[L - 1][1];++col){

         if(sparse[L - 1][row*dim[L - 1][1] + col]){

            cout << endl;
            cout << row << "\t" << col << endl;
            cout << endl;
            cout << HamOp::ostates[L - 2][row] << endl;
            cout << HamOp::ostates[0][0] << endl;
            cout << endl;
            cout << *op[L - 1][row*dim[L - 1][1] + col] << endl;
            cout << endl;

         }

      }


}

/**
 * access to the operators from outside of the class
 */
const QSDArray<2> &Operator::gop(int site,int row,int col){

   return *op[site][row*dim[site][1] + col];

}

/**
 * access to the sparsity info from outside of the class
 */
bool Operator::gsparse(int site,int row,int col){

   return sparse[site][row*dim[site][1] + col];

}
