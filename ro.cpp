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
    * construct a renormalized operator for <wccd| H T |wccd> dmrg style, write everything to disk
    * @param dir left renormalized or right renormalized
    * @param wccd input MPS
    */
   void construct(const MPS_DIRECTION &dir,const MPO<Quantum> &T,const MPS<Quantum> &wccd){

      int L = wccd.size();

      if(dir == mpsxx::Left){

         QSDArray<5> tmp5;

         enum {j,k,l,m,n,o};

         QSDindexed_contract(1.0,wccd[0],shape(j,k,l),T[0],shape(m,n,k,o),0.0,tmp5,shape(j,m,n,l,o));

         //merge the row quantumnumbers together
         TVector<Qshapes<Quantum>,2> qmerge;
         TVector<Dshapes,2> dmerge;
         
         for(int i = 0;i < 2;++i){

            qmerge[i] = tmp5.qshape(i);
            dmerge[i] = tmp5.dshape(i);

         }

         QSTmergeInfo<2> info(qmerge,dmerge);

         //then merge into tmp4
         QSDArray<4> tmp4;
         QSTmerge(info,tmp5,tmp4);

         //first site:
         for(int col = 0;col < Operator::gdim(0,1);++col)
            print_op(mpsxx::Left,0,col,Operator::gop(0,0,col),tmp4,wccd[0]);

         //middle sites:
         for(int i = 1;i < L - 1;++i){

            cout << endl;
            cout << i << endl;
            cout << endl;

            //get the previous operator from memory and paste the next sites onto it:
            tmp5.clear();

            //first row = 0: id coming in
            get_op(mpsxx::Left,i - 1,0,T[i],wccd[i],tmp5);

            int col = 0;

            while(HamOp::ostates[i][col].size() == 1){

               if(Operator::gsparse(i,0,col))
                  print_op(mpsxx::Left,i,col,Operator::gop(i,0,col),tmp5);

               ++col;

            }

            while(HamOp::ostates[i][col].size() == 2){

               if(Operator::gsparse(i,0,col))
                  print_op(mpsxx::Left,i,col,Operator::gop(i,0,col),tmp5);

               ++col;

            }

            int osbar = col;

            std::vector< QSDArray<3> > os(HamOp::ostates[i].size() - osbar);

            while(col < HamOp::ostates[i].size()){

               if(Operator::gsparse(i,0,col))
                  QSDcontract(1.0,tmp5,shape(3,4),Operator::gop(i,0,col),shape(0,1),0.0,os[col - osbar]);

               ++col;

            }

            for(int row = 1;row < HamOp::ostates[i - 1].size();++row){

               cout << i << "\t" << row << "\t" << HamOp::ostates[i-1].size() << endl;

               tmp5.clear();
               get_op(mpsxx::Left,i - 1,row,T[i],wccd[i],tmp5);

               for(int col = 0;col < osbar;++col)
                  if(Operator::gsparse(i,row,col))
                     print_op(mpsxx::Left,i,col,Operator::gop(i,row,col),tmp5);

               for(int col = osbar;col < HamOp::ostates[i].size();++col)
                  if(Operator::gsparse(i,row,col))
                    QSDcontract(1.0,tmp5,shape(3,4),Operator::gop(i,row,col),shape(0,1),1.0,os[col - osbar]);

            }

            //print out the os part for the next site
            for(int col = osbar;col < HamOp::ostates[i].size();++col)
               save(mpsxx::Left,i,col,os[col - osbar]);

         }

      }
      else{

         QSDArray<5> tmp5;

         enum {j,k,l,m,n,o};

         QSDindexed_contract(1.0,wccd[L - 1],shape(j,k,l),T[L - 1],shape(m,n,k,o),0.0,tmp5,shape(j,m,n,l,o));

         //merge the column quantumnumbers together
         TVector<Qshapes<Quantum>,2> qmerge;
         TVector<Dshapes,2> dmerge;

         for(int i = 0;i < 2;++i){

            qmerge[i] = tmp5.qshape(3 + i);
            dmerge[i] = tmp5.dshape(3 + i);

         }

         QSTmergeInfo<2> info(qmerge,dmerge);

         //then merge into tmp4
         QSDArray<4> tmp4;
         QSTmerge(tmp5,info,tmp4);

         std::vector< std::vector<bool> > flag(2);

         flag[0].resize(HamOp::ostates[L - 2].size());

         //last site:
         for(int row = 0;row < HamOp::ostates[L - 2].size();++row){

            flag[0][row] = false;

            if(Operator::gsparse(L-1,row,0)){

               print_op(mpsxx::Right,L-1,row,Operator::gop(L-1,row,0),tmp4,wccd[L-1]);
               flag[0][row] = true;

            }

         }

         //rest of the sites
         for(int i = L - 2;i > 0;--i){

            //copy the flag and init new one
            flag[1] = flag[0];

            flag[0].resize(HamOp::ostates[i - 1].size());

            for(int row = 0;row < flag[0].size();++row)
               flag[0][row] = false;

            cout << endl;
            cout << i << endl;
            cout << endl;

            //first find where the incoming singles stop
            int isbar = 0;

            while(HamOp::ostates[i - 1][isbar].size() == 1)
               isbar++;

            std::vector< QSDArray<3> > is(isbar);

            int col = 0;

            while(HamOp::ostates[i][col].size() == 1){

               cout << i << "\t" << col << "\t" << HamOp::ostates[i].size() << endl;

               tmp5.clear();
               get_op(mpsxx::Right,i + 1,col,T[i],wccd[i],tmp5);

               for(int row = 0;row < isbar;++row)
                  if(Operator::gsparse(i,row,col))
                     QSDcontract(1.0,tmp5,shape(3,4),Operator::gop(i,row,col),shape(0,1),1.0,is[row]);

               ++col;

            }

            while(HamOp::ostates[i][col].size() == 2){

               cout << i << "\t" << col << "\t" << HamOp::ostates[i].size() << endl;

               tmp5.clear();

               if(flag[1][col]){

                  get_op(mpsxx::Right,i + 1,col,T[i],wccd[i],tmp5);

                  for(int row = 0;row < isbar;++row)
                     if(Operator::gsparse(i,row,col))
                        QSDcontract(1.0,tmp5,shape(3,4),Operator::gop(i,row,col),shape(0,1),1.0,is[row]);

                  for(int row = isbar;row < HamOp::ostates[i - 1].size();++row){

                     if(Operator::gsparse(i,row,col)){

                        print_op(mpsxx::Right,i,row,Operator::gop(i,row,col),tmp5);
                        flag[0][row] = true;

                     }

                  }

               }

               ++col;

            }

            int osbar = col;

            //read in the last columns : closed and outgoing singles coming in
            std::vector< QSDArray<5> > os(HamOp::ostates[i].size() - osbar);

            for(int col = osbar;col < HamOp::ostates[i].size();++col)
               get_op(mpsxx::Right,i + 1,col,T[i],wccd[i],os[col - osbar]);

            for(int row = 0;row < isbar;++row){

               cout << i << "\t" << row << "\t" << HamOp::ostates[i - 1].size() << endl;

               for(int col = osbar;col < HamOp::ostates[i].size();++col)
                  if(Operator::gsparse(i,row,col))
                     QSDcontract(1.0,os[col - osbar],shape(3,4),Operator::gop(i,row,col),shape(0,1),1.0,is[row]);

               save(mpsxx::Right,i,row,is[row]);
               flag[0][row] = true;

            }

            //doubles
            int row = isbar;

            QSDArray<3> tmp3;

            while(HamOp::ostates[i - 1][row].size() == 2){

               cout << i << "\t" << row << "\t" << HamOp::ostates[i - 1].size() << endl;

               tmp3.clear();

               if(flag[0][row])
                  read(mpsxx::Right,i,row,tmp3);

               for(int col = osbar;col < HamOp::ostates[i].size();++col)
                  if(Operator::gsparse(i,row,col))
                     QSDcontract(1.0,os[col - osbar],shape(3,4),Operator::gop(i,row,col),shape(0,1),1.0,tmp3);

               save(mpsxx::Right,i,row,tmp3);
               flag[0][row] = true;

               ++row;

            }

            while(row < HamOp::ostates[i - 1].size()){

               cout << i << "\t" << row << "\t" << HamOp::ostates[i - 1].size() << endl;

               for(int col = osbar;col < HamOp::ostates[i].size();++col)
                  if(Operator::gsparse(i,row,col))
                     print_op(mpsxx::Right,i,row,Operator::gop(i,row,col),os[col - osbar]);

               ++row;

            }

         }
         
      }

   }

   /**
    * print the contraction of two QSDArrays A on site 'site'
    */
   void print_op(const MPS_DIRECTION &dir,int site,int opnum,const QSDArray<2> &op,const QSDArray<4> &TA,const QSDArray<3> &A){

      if(dir == mpsxx::Left){

         enum {j,k,l,m,n};

         QSDArray<3> tmp;

         QSDindexed_contract(1.0,A.conjugate(),shape(j,k,m),op,shape(k,n),0.0,tmp,shape(j,n,m));

         QSDArray<3> E;

         QSDindexed_contract(1.0,TA,shape(j,n,l,k),tmp,shape(j,n,m),0.0,E,shape(l,k,m));

         char name[100];

         sprintf(name,"scratch/Left/site_%d/%d.mpx",site,opnum);

         std::ofstream fout(name);
         boost::archive::binary_oarchive oar(fout);

         oar << E;

      }
      else{

         enum {j,k,l,m,n};

         QSDArray<3> tmp;

         QSDindexed_contract(1.0,A.conjugate(),shape(m,k,j),op,shape(k,n),0.0,tmp,shape(m,n,j));

         QSDArray<3> E;

         QSDindexed_contract(1.0,TA,shape(l,k,n,j),tmp,shape(m,n,j),0.0,E,shape(l,k,m));

         char name[100];

         sprintf(name,"scratch/Right/site_%d/%d.mpx",site,opnum);

         std::ofstream fout(name);
         boost::archive::binary_oarchive oar(fout);

         oar << E;

      }

   }

   /**
    * read in the operator on the previous site with number opnum, paste the next site to it:
    */
   void get_op(const MPS_DIRECTION &dir,int site,int opnum,const QSDArray<4> &T,const QSDArray<3> &A,QSDArray<5> &E_op){

      if(dir == mpsxx::Left){

         enum {i,j,k,l,m,n,o};

         char name[100];

         sprintf(name,"scratch/Left/site_%d/%d.mpx",site,opnum);

         std::ifstream fin(name);
         boost::archive::binary_iarchive iar(fin);

         QSDArray<3> tmp3;

         iar >> tmp3;

         QSDArray<4> tmp4;

         QSDindexed_contract(1.0,tmp3,shape(o,k,m),A,shape(o,n,l),0.0,tmp4,shape(l,k,m,n));

         QSDArray<4> tmp4_bis;

         QSDindexed_contract(1.0,tmp4,shape(l,o,m,n),T,shape(o,j,n,k),0.0,tmp4_bis,shape(l,k,m,j));

         QSDindexed_contract(1.0,tmp4_bis,shape(l,k,o,j),A.conjugate(),shape(o,i,m),0.0,E_op,shape(l,k,m,i,j));

      }
      else{

         enum {i,j,k,l,m,n,o};

         char name[100];

         sprintf(name,"scratch/Right/site_%d/%d.mpx",site,opnum);

         std::ifstream fin(name);
         boost::archive::binary_iarchive iar(fin);

         QSDArray<3> tmp3;

         iar >> tmp3;

         QSDArray<4> tmp4;

         QSDindexed_contract(1.0,A,shape(o,n,l),tmp3,shape(l,k,m),0.0,tmp4,shape(o,k,m,n));

         QSDArray<4> tmp4_bis;

         QSDindexed_contract(1.0,tmp4,shape(o,k,m,n),T,shape(l,j,n,k),0.0,tmp4_bis,shape(o,l,m,j));

         QSDindexed_contract(1.0,A.conjugate(),shape(m,i,n),tmp4_bis,shape(l,k,n,j),0.0,E_op,shape(l,k,m,i,j));

      }

   }

   /**
    * print the contraction of two QSDArrays A on site 'site'
    */
   void print_op(const MPS_DIRECTION &dir,int site,int opnum,const QSDArray<2> &op,const QSDArray<5> &E_op){

      enum {i,j,k,l,m};

      QSDArray<3> E;
      E.clear();

      QSDcontract(1.0,E_op,shape(3,4),op,shape(0,1),0.0,E);

      char name[100];

      if(dir == mpsxx::Left)
         sprintf(name,"scratch/Left/site_%d/%d.mpx",site,opnum);
      else
         sprintf(name,"scratch/Right/site_%d/%d.mpx",site,opnum);

      std::ofstream fout(name);
      boost::archive::binary_oarchive oar(fout);

      oar << E;

   }

   /**
    * read in the matrix on site i with identifier opnum
    */
   void read(const MPS_DIRECTION &dir,int site,int opnum,QSDArray<3> &E){

      char name[100];

      if(dir == mpsxx::Left)
         sprintf(name,"scratch/Left/site_%d/%d.mpx",site,opnum);
      else
         sprintf(name,"scratch/Right/site_%d/%d.mpx",site,opnum);

      std::ifstream fin(name);
      boost::archive::binary_iarchive iar(fin);

      iar >> E;

   }


   /**
    * save the matrix on site i with id opnum
    */
   void save(const MPS_DIRECTION &dir,int site,int opnum,const QSDArray<3> &E){

      char name[100];

      if(dir == mpsxx::Left)
         sprintf(name,"scratch/Left/site_%d/%d.mpx",site,opnum);
      else
         sprintf(name,"scratch/Right/site_%d/%d.mpx",site,opnum);

      std::ofstream fout(name);
      boost::archive::binary_oarchive oar(fout);

      oar << E;

   }

}
