#include <iostream>
#include <fstream>
#include <cmath>

using std::endl;
using std::ostream;
using std::ofstream;
using std::ifstream;

#include "include.h"

using namespace btas;
using namespace mpsxx;

namespace ro {

   /**
    * construct a renormalized operator for <B|O|A> and store it in an MPS
    * @param dir left renormalized or right renormalized
    * @param A input MPS
    * @param O input MPO
    * @param B input MPS
    */
   MPS<Quantum> construct(const MPS_DIRECTION &dir,const MPS<Quantum> &A,const MPO<Quantum> &O,const MPS<Quantum> &B){

      int L = O.size();

      MPS<Quantum> RO(L - 1);

      if(dir == mpsxx::Left){

         enum {j,k,l,m,n,o};

         //from left to right
         QSDArray<5> loc;

         QSDindexed_contract(1.0,O[0],shape(m,n,k,o),A[0],shape(j,k,l),0.0,loc,shape(j,m,n,l,o));

         //merge 2 rows together
         TVector<Qshapes<Quantum>,2> qmerge;
         TVector<Dshapes,2> dmerge;

         for(int i = 0;i < 2;++i){

            qmerge[i] = loc.qshape(i);
            dmerge[i] = loc.dshape(i);

         }

         QSTmergeInfo<2> info(qmerge,dmerge);

         QSDArray<4> tmp;
         QSTmerge(info,loc,tmp);

         //this will contain the right going part
         QSDindexed_contract(1.0,B[0].conjugate(),shape(j,k,l),tmp,shape(j,k,m,n),0.0,RO[0],shape(m,n,l));

         QSDArray<4> I1;
         QSDArray<4> I2;

         for(int i = 1;i < L - 1;++i){

            I1.clear();

            QSDindexed_contract(1.0,RO[i - 1],shape(j,k,l),A[i],shape(j,m,n),0.0,I1,shape(k,l,n,m));

            I2.clear();

            QSDindexed_contract(1.0,I1,shape(k,l,n,m),O[i],shape(k,j,m,o),0.0,I2,shape(l,j,n,o));

            QSDindexed_contract(1.0,I2,shape(l,j,n,o),B[i].conjugate(),shape(l,j,k),0.0,RO[i],shape(n,o,k));

         }

      }
      else{

         enum {j,k,l,m,n,o};

         //from right to left
         QSDArray<5> loc;

         QSDindexed_contract(1.0,O[L - 1],shape(j,k,l,m),A[L - 1],shape(o,l,n),0.0,loc,shape(o,j,k,n,m));

         //merge 2 columns together
         TVector<Qshapes<Quantum>,2> qmerge;
         TVector<Dshapes,2> dmerge;

         for(int i = 0;i < 2;++i){

            qmerge[i] = loc.qshape(3 + i);
            dmerge[i] = loc.dshape(3 + i);

         }

         QSTmergeInfo<2> info(qmerge,dmerge);

         QSDArray<4> tmp;
         QSTmerge(loc,info,tmp);

         //this will contain the left going part
         QSDindexed_contract(1.0,tmp,shape(j,k,l,m),B[L-1].conjugate(),shape(n,l,m),0.0,RO[L - 2],shape(j,k,n));

         QSDArray<4> I1;
         QSDArray<4> I2;

         for(int i = L - 2;i > 0;--i){

            I1.clear();

            QSDindexed_contract(1.0,A[i],shape(j,k,l),RO[i],shape(l,m,n),0.0,I1,shape(j,k,m,n));

            I2.clear();

            QSDindexed_contract(1.0,O[i],shape(l,o,k,m),I1,shape(j,k,m,n),0.0,I2,shape(j,l,o,n));

            QSDindexed_contract(1.0,B[i].conjugate(),shape(k,o,n),I2,shape(j,l,o,n),0.0,RO[i - 1],shape(j,l,k));

         }

      }

      return RO;

   }

   /**
    * check if the renormalized operators are correctly constructed
    */
   void check(const MPS<Quantum> &ror,const MPS<Quantum> &rol){

      for(int i = 0;i < ror.size();++i)
         cout << QSDdotc(rol[i],ror[i].conjugate()) << endl;

   }

   /**
    * check if the renormalized operators are correctly constructed
    */
   void check(const MPO<Quantum> &ror,const MPO<Quantum> &rol){

      for(int i = 0;i < ror.size();++i)
         cout << QSDdotc(rol[i],ror[i].conjugate()) << endl;

   }


   /**
    * construct a renormalized operator for <tccd| T H |wccd> and store it in an MPO
    * @param dir left renormalized or right renormalized
    * @param tccd input MPS
    * @param T input MPO
    * @param H input MPO
    * @param wccd input MPS
    */
   MPO<Quantum> construct(const MPS_DIRECTION &dir,const MPS<Quantum> &tccd,const MPO<Quantum> &T,const MPO<Quantum> &H,const MPS<Quantum> &wccd){

      int L = wccd.size();

      MPO<Quantum> RO(L - 1);

      if(dir == mpsxx::Left){

         enum {j,k,l,m,n,o,p};

         //from left to right
         QSDArray<5> tmp5;

         QSDindexed_contract(1.0,T[0],shape(m,n,k,o),tccd[0],shape(j,k,l),0.0,tmp5,shape(j,m,n,l,o));

         //merge 2 rows together
         TVector<Qshapes<Quantum>,2> qmerge;
         TVector<Dshapes,2> dmerge;

         for(int i = 0;i < 2;++i){

            qmerge[i] = tmp5.qshape(i);
            dmerge[i] = tmp5.dshape(i);

         }

         QSTmergeInfo<2> info(qmerge,dmerge);

         QSDArray<4> tmp4;
         QSTmerge(info,tmp5,tmp4);

         QSDArray<6> tmp6;

         QSDindexed_contract(1.0,H[0],shape(j,m,n,p),tmp4,shape(k,n,l,o),0.0,tmp6,shape(k,j,m,l,o,p));

         //merge again
         for(int i = 0;i < 2;++i){

            qmerge[i] = tmp6.qshape(i);
            dmerge[i] = tmp6.dshape(i);

         }

         info.reset(qmerge,dmerge);

         tmp5.clear();
         QSTmerge(info,tmp6,tmp5);

         //this will contain the right going part
         QSDindexed_contract(1.0,wccd[0].conjugate(),shape(j,m,k),tmp5,shape(j,m,l,o,p),0.0,RO[0],shape(l,o,p,k));

         QSDArray<5> I1;
         QSDArray<5> I2;

         for(int i = 1;i < L - 1;++i){

            tmp5.clear();

            QSDindexed_contract(1.0,RO[i - 1],shape(j,k,l,m),tccd[i],shape(j,n,o),0.0,tmp5,shape(o,n,k,l,m));

            I1.clear();

            QSDindexed_contract(1.0,tmp5,shape(o,n,k,l,m),T[i],shape(k,j,n,p),0.0,I1,shape(o,p,j,l,m));

            I2.clear();

            QSDindexed_contract(1.0,I1,shape(o,p,j,l,m),H[i],shape(l,k,j,n),0.0,I2,shape(o,p,n,k,m));

            QSDindexed_contract(1.0,I2,shape(o,p,n,k,m),wccd[i].conjugate(),shape(m,k,l),0.0,RO[i],shape(o,p,n,l));

         }

      }
      else{

         enum {j,k,l,m,n,o,p};

         //from right to left
         QSDArray<5> tmp5;

         QSDindexed_contract(1.0,T[L - 1],shape(m,n,k,o),tccd[L - 1],shape(j,k,l),0.0,tmp5,shape(j,m,n,l,o));

         //merge 2 columns together
         TVector<Qshapes<Quantum>,2> qmerge;
         TVector<Dshapes,2> dmerge;

         for(int i = 0;i < 2;++i){

            qmerge[i] = tmp5.qshape(i + 3);
            dmerge[i] = tmp5.dshape(i + 3);

         }

         QSTmergeInfo<2> info(qmerge,dmerge);

         QSDArray<4> tmp4;
         QSTmerge(tmp5,info,tmp4);

         QSDArray<6> tmp6;

         QSDindexed_contract(1.0,H[L - 1],shape(k,l,n,o),tmp4,shape(j,m,n,p),0.0,tmp6,shape(j,m,k,l,p,o));

         //merge again
         for(int i = 0;i < 2;++i){

            qmerge[i] = tmp6.qshape(i + 4);
            dmerge[i] = tmp6.dshape(i + 4);

         }

         info.reset(qmerge,dmerge);

         tmp5.clear();
         QSTmerge(tmp6,info,tmp5);

         //this will contain the right going part
         QSDindexed_contract(1.0,wccd[L - 1].conjugate(),shape(n,l,p),tmp5,shape(j,m,k,l,p),0.0,RO[L - 2],shape(j,m,k,n));

         QSDArray<5> I1;
         QSDArray<5> I2;

         for(int i = L - 2;i > 0;--i){

            tmp5.clear();

            QSDindexed_contract(1.0,tccd[i],shape(n,o,j),RO[i],shape(j,k,l,m),0.0,tmp5,shape(n,o,k,l,m));

            I1.clear();

            QSDindexed_contract(1.0,T[i],shape(j,p,o,k),tmp5,shape(n,o,k,l,m),0.0,I1,shape(n,j,p,l,m));

            I2.clear();

            QSDindexed_contract(1.0,H[i],shape(k,o,p,l),I1,shape(n,j,p,l,m),0.0,I2,shape(n,j,k,o,m));

            QSDindexed_contract(1.0,wccd[i].conjugate(),shape(l,o,m),I2,shape(n,j,k,o,m),0.0,RO[i - 1],shape(n,j,k,l));

         }

      }

      return RO;

   }

   /**
    * construct a renormalized operator for <wccd| H |wccd> dmrg style, write everything to disk
    * @param dir left renormalized or right renormalized
    * @param wccd input MPS
    */
   void construct(const MPS_DIRECTION &dir,const MPS<Quantum> &wccd){

      int L = wccd.size();
     
      //print id
      print_op(0,0,Operator::id,wccd[0]);

      //insert singles
      print_op(0,1,Operator::cus,wccd[0]);
      print_op(0,2,Operator::cd,wccd[0]);
      print_op(0,3,Operator::aus,wccd[0]);
      print_op(0,4,Operator::ad,wccd[0]);

      //insert doubles
      print_op(0,5,Operator::cucd,wccd[0]);
      print_op(0,6,Operator::cuau,wccd[0]);
      print_op(0,7,Operator::cuad,wccd[0]);
      print_op(0,8,Operator::cdau,wccd[0]);
      print_op(0,9,Operator::cdad,wccd[0]);
      print_op(0,10,Operator::auad,wccd[0]);//wrong sign

      int col = 11;

      //insert triples: construct complementary operator
      for(int j = 1;j < L;++j){

         print_op(0,col,Operator::tcuf[j-1],wccd[0]);col++;
         print_op(0,col,Operator::tcdf[j-1],wccd[0]);col++;
         print_op(0,col,Operator::tauf[j-1],wccd[0]);col++;
         print_op(0,col,Operator::tadf[j-1],wccd[0]);col++;

      }

      //last term, local term already closed:
      print_op(0,col,Operator::local,wccd[0]);

      //middle sites
      for(int i = 1;i < 2;++i){

/*
         //insert singles
         insert_crea_up_s(mpo[i],0,1,1.0);
         insert_crea_down(mpo[i],0,2,1.0);
         insert_anni_up_s(mpo[i],0,3,1.0);
         insert_anni_down(mpo[i],0,4,1.0);

         //insert signs
         row = 1;
         column = 5;

         while(HamOp::ostates[i-1][row].size() == 1){

            insert_sign(mpo[i],row,column);
            ++row;
            ++column;

         }

         //insert doubles
         insert_crea_up_crea_down(mpo[i],0,column,1.0);++column;
         insert_crea_up_anni_up(mpo[i],0,column,1.0);++column;
         insert_crea_up_anni_down(mpo[i],0,column,1.0);++column;
         insert_crea_down_anni_up(mpo[i],0,column,1.0);++column;
         insert_crea_down_anni_down(mpo[i],0,column,1.0);++column;
         insert_anni_down_anni_up(mpo[i],0,column,1.0);++column;

         //insert singles to form outgoing doubles
         row = 1;

         while(HamOp::ostates[i-1][row].size() == 1){//create up

            insert_crea_up(mpo[i],row,column,1.0);
            ++row;
            ++column;

         }

         row = 1;

         while(HamOp::ostates[i-1][row].size() == 1){//create down with sign

            insert_crea_down_s(mpo[i],row,column,1.0);
            ++row;
            ++column;

         }

         row = 1;

         while(HamOp::ostates[i-1][row].size() == 1){//annihilate up

            insert_anni_up(mpo[i],row,column,1.0);
            ++row;
            ++column;

         }

         row = 1;

         while(HamOp::ostates[i-1][row].size() == 1){//annihilate down with sign

            insert_anni_down_s(mpo[i],row,column,1.0);
            ++row;
            ++column;

         }

         //copy the doubles from previous tensor: identity
         while(HamOp::ostates[i-1][row].size() == 2){

            insert_id(mpo[i],row,column);
            ++row;
            ++column;

         }

         //HERE STARTS THE COMPLEMENTARY OPERATOR STUFF!
         while(column < HamOp::ostates[i].size() - 1){

            int j = HamOp::ostates[i][column].gsite(0);
            int sj = HamOp::ostates[i][column].gspin(0);
            int aj = HamOp::ostates[i][column].gact(0);

            //first row
            if(sj == 0 && aj == 0)
               insert_triple_crea_up_first(mpo[i],0,column,t(i,j),V(i,i,j,i));
            else if(sj == 1 && aj == 0)
               insert_triple_crea_down_first(mpo[i],0,column,t(i,j),V(i,i,j,i));
            else if(sj == 0 && aj == 1)
               insert_triple_anni_up_first(mpo[i],0,column,t(i,j),V(i,i,j,i));
            else
               insert_triple_anni_down_first(mpo[i],0,column,t(i,j),V(i,i,j,i));

            //singles coming in:
            row = 1;

            while(HamOp::ostates[i-1][row].size() == 1){

               std::vector<double> val(2);

               std::vector<int> v = Ostate::get_double_complement(i,HamOp::ostates[i-1][row],HamOp::ostates[i][column],V,val);

               if(v.size() == 1){

                  if(v[0] == 0)
                     insert_anni_down_anni_up(mpo[i],row,column,-val[0]);//extra minus sign!
                  else if(v[0] == 1)
                     insert_crea_up_crea_down(mpo[i],row,column,val[0]);
                  else if(v[0] == 2)
                     insert_crea_up_anni_down(mpo[i],row,column,val[0]);
                  else
                     insert_crea_down_anni_up(mpo[i],row,column,val[0]);

               }
               else if(v.size() == 2)//with sign because in the middle!
                  insert_pair_s(mpo[i],row,column,val);

               ++row;

            }

            //pairs coming in
            while(HamOp::ostates[i-1][row].size() == 2){

               double val;

               std::vector<int> v = Ostate::get_single_complement(i,HamOp::ostates[i-1][row],HamOp::ostates[i][column],V,val);

               if(v.size() > 0){

                  if(v[0] == 0 && v[1] == 0)//create spin up
                     insert_crea_up_s(mpo[i],row,column,val);
                  else if(v[0] == 1 && v[1] == 0)//create spin down
                     insert_crea_down(mpo[i],row,column,val);
                  else if(v[0] == 0 && v[1] == 1)//annihilate spin up
                     insert_anni_up_s(mpo[i],row,column,val);
                  else//annihilate spin down
                     insert_anni_down(mpo[i],row,column,val);

               }

               ++row;

            }

            //signs: find row and column which are connected
            while(HamOp::ostates[i-1][row] != HamOp::ostates[i][column])
               ++row;

            insert_sign(mpo[i],row,column);

            column++;

         }

         //last column! closing down everything!

         //first row
         insert_local(mpo[i],0,column,t(i,i),V(i,i,i,i));

         //close down the singles coming in with a triplet
         row = 1;

         while(HamOp::ostates[i-1][row].size() == 1){

            //incoming operator
            int k = HamOp::ostates[i-1][row].gsite(0);
            int sk = HamOp::ostates[i-1][row].gspin(0);
            int ak = HamOp::ostates[i-1][row].gact(0);

            if(sk == 0 && ak == 0)//create up coming in
               insert_triple_crea_up_last(mpo[i],row,column,V(k,i,i,i));
            else if(sk == 1 && ak == 0)//create down coming in
               insert_triple_crea_down_last(mpo[i],row,column,-V(i,k,i,i));
            else if(sk == 0 && ak == 1)//annihilate up coming in
               insert_triple_anni_up_last(mpo[i],row,column,V(i,i,k,i));
            else //annihilate down coming in
               insert_triple_anni_down_last(mpo[i],row,column,-V(i,i,i,k));

            ++row;

         }

         //close down the doubles coming in with a pair
         while(HamOp::ostates[i-1][row].size() == 2){

            std::vector<double> val(2);

            std::vector<int> v = Ostate::get_closing_pair_in(i,HamOp::ostates[i-1][row],V,val);

            if(v.size() == 1){

               if(v[0] == 0)
                  insert_anni_down_anni_up(mpo[i],row,column,-val[0]);//extra minus sign!
               else if(v[0] == 1)
                  insert_crea_up_crea_down(mpo[i],row,column,val[0]);
               else if(v[0] == 2)
                  insert_crea_up_anni_down(mpo[i],row,column,val[0]);
               else
                  insert_crea_down_anni_up(mpo[i],row,column,val[0]);

            }
            else if(v.size() == 2)
               insert_pair(mpo[i],row,column,val);

            ++row;

         }

         //close down the complementary operators of this site
         //basically the first 4 incoming states should be closed down
         insert_crea_up(mpo[i],row,column,1.0);
         ++row;

         insert_crea_down_s(mpo[i],row,column,1.0);
         ++row;

         insert_anni_up(mpo[i],row,column,1.0);
         ++row;

         insert_anni_down_s(mpo[i],row,column,1.0);
         ++row;

         //finally insert the identity on the lower right element
         insert_id(mpo[i],HamOp::ostates[i-1].size() - 1,HamOp::ostates[i].size() - 1);
*/
      }


   }

   /**
    * print the contraction of two QSDArrays A on site 'site'
    */
   void print_op(int site,int opnum,const QSDArray<2> &op,const QSDArray<3> &A){

      if(site == 0){

         enum {j,k,l,m};

         QSDArray<3> tmp;

         QSDindexed_contract(1.0,A,shape(j,k,l),op,shape(m,k),0.0,tmp,shape(j,m,l));

         QSDArray<2> E;

         QSDindexed_contract(1.0,tmp,shape(j,k,l),A.conjugate(),shape(j,k,m),0.0,E,shape(l,m));

         cout << E << endl;

         char name[100];

         sprintf(name,"scratch/site_%d/%d.mpx",site,opnum);

         std::ofstream fout(name);
         boost::archive::binary_oarchive oar(fout);

         oar << E;

      }

   }

}
