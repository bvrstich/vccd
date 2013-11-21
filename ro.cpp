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

      if(dir == mpsxx::Left){

         //first site:
         for(int col = 0;col < Operator::gdim(0,1);++col)
            print_op(mpsxx::Left,0,col,Operator::gop(0,0,col),wccd[0]);

         QSDArray<2> tmp2;

         //middle sites:
         for(int i = 1;i < L - 1;++i){

            cout << endl;
            cout << i << endl;
            cout << endl;

            //get the previous operator from memory and paste the next sites onto it:
            QSDArray<4> tmp4;

            //first row = 0: id coming in
            get_op(mpsxx::Left,i - 1,0,wccd[i],tmp4);

            int col = 0;

            while(HamOp::ostates[i][col].size() == 1){

               if(Operator::gsparse(i,0,col))
                  print_op(mpsxx::Left,i,col,Operator::gop(i,0,col),tmp4);

               ++col;

            }

            while(HamOp::ostates[i][col].size() == 2){

               if(Operator::gsparse(i,0,col))
                  print_op(mpsxx::Left,i,col,Operator::gop(i,0,col),tmp4);

               ++col;

            }

            int osbar = col;

            std::vector< QSDArray<2> > os(HamOp::ostates[i].size() - osbar);

            while(col < HamOp::ostates[i].size()){

               if(Operator::gsparse(i,0,col))
                  QSDcontract(1.0,tmp4,shape(3,4),Operator::gop(i,0,col),shape(0,1),0.0,os[col - osbar]);

               ++col;

            }

            for(int row = 1;row < HamOp::ostates[i - 1].size();++row){

               cout << i << "\t" << row << "\t" << HamOp::ostates[i-1].size() << endl;

               tmp4.clear();
               get_op(mpsxx::Left,i - 1,row,wccd[i],tmp4);

               for(int col = 0;col < osbar;++col)
                  if(Operator::gsparse(i,row,col))
                     print_op(mpsxx::Left,i,col,Operator::gop(i,row,col),tmp4);

               for(int col = osbar;col < HamOp::ostates[i].size();++col)
                  if(Operator::gsparse(i,row,col))
                     QSDcontract(1.0,tmp4,shape(3,4),Operator::gop(i,row,col),shape(0,1),1.0,os[col - osbar]);

            }

            //print out the os part for the next site
            for(int col = osbar;col < HamOp::ostates[i].size();++col)
               save(mpsxx::Left,i,col,os[col - osbar]);

         }

      }
      else{

         std::vector< std::vector<bool> > flag(2);

         flag[0].resize(HamOp::ostates[L - 2].size());

         //last site:
         for(int row = 0;row < HamOp::ostates[L - 2].size();++row){

            flag[0][row] = false;

            if(Operator::gsparse(L-1,row,0)){

               print_op(mpsxx::Right,L-1,row,Operator::gop(L-1,row,0),wccd[L-1]);
               flag[0][row] = true;

            }

         }

         QSDArray<2> tmp2;

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

            std::vector< QSDArray<2> > is(isbar);

            QSDArray<4> tmp4;

            int col = 0;

            while(HamOp::ostates[i][col].size() == 1){

               cout << i << "\t" << col << "\t" << HamOp::ostates[i].size() << endl;

               tmp4.clear();
               get_op(mpsxx::Right,i + 1,col,wccd[i],tmp4);

               for(int row = 0;row < isbar;++row)
                  if(Operator::gsparse(i,row,col))
                     QSDcontract(1.0,tmp4,shape(3,4),Operator::gop(i,row,col),shape(0,1),1.0,is[row]);

               ++col;

            }

            while(HamOp::ostates[i][col].size() == 2){

               cout << i << "\t" << col << "\t" << HamOp::ostates[i].size() << endl;

               tmp4.clear();

               if(flag[1][col]){

                  get_op(mpsxx::Right,i + 1,col,wccd[i],tmp4);

                  for(int row = 0;row < isbar;++row)
                     if(Operator::gsparse(i,row,col))
                        QSDcontract(1.0,tmp4,shape(3,4),Operator::gop(i,row,col),shape(0,1),1.0,is[row]);

                  for(int row = isbar;row < HamOp::ostates[i - 1].size();++row){

                     if(Operator::gsparse(i,row,col)){

                        print_op(mpsxx::Right,i,row,Operator::gop(i,row,col),tmp4);
                        flag[0][row] = true;

                     }

                  }

               }

               ++col;

            }

            int osbar = col;

            //read in the last columns : closed and outgoing singles coming in
            std::vector< QSDArray<4> > os(HamOp::ostates[i].size() - osbar);

            for(int col = osbar;col < HamOp::ostates[i].size();++col)
               get_op(mpsxx::Right,i + 1,col,wccd[i],os[col - osbar]);

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

            while(HamOp::ostates[i - 1][row].size() == 2){

               cout << i << "\t" << row << "\t" << HamOp::ostates[i - 1].size() << endl;

               tmp2.clear();

               if(flag[0][row])
                  read(mpsxx::Right,i,row,tmp2);

               for(int col = osbar;col < HamOp::ostates[i].size();++col)
                  if(Operator::gsparse(i,row,col))
                     QSDcontract(1.0,os[col - osbar],shape(3,4),Operator::gop(i,row,col),shape(0,1),1.0,tmp2);

               save(mpsxx::Right,i,row,tmp2);
               flag[0][row] = true;

               ++row;

            }

            while(row < HamOp::ostates[i - 1].size()){

               cout << i << "\t" << row << "\t" << HamOp::ostates[i - 1].size() << endl;

               tmp2.clear();

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
   void print_op(const MPS_DIRECTION &dir,int site,int opnum,const QSDArray<2> &op,const QSDArray<3> &A){

      if(dir == mpsxx::Left){

         enum {j,k,l,m};

         QSDArray<3> tmp;

         QSDindexed_contract(1.0,A,shape(j,k,l),op,shape(m,k),0.0,tmp,shape(j,m,l));

         QSDArray<2> E;

         QSDindexed_contract(1.0,tmp,shape(j,k,l),A.conjugate(),shape(j,k,m),0.0,E,shape(l,m));

         char name[100];

         sprintf(name,"scratch/Left/site_%d/%d.mpx",site,opnum);

         std::ofstream fout(name);
         boost::archive::binary_oarchive oar(fout);

         oar << E;

      }
      else{

         enum {j,k,l,m};

         QSDArray<3> tmp;

         QSDindexed_contract(1.0,A,shape(l,k,j),op,shape(m,k),0.0,tmp,shape(l,m,j));

         QSDArray<2> E;

         QSDindexed_contract(1.0,tmp,shape(l,k,j),A.conjugate(),shape(m,k,j),0.0,E,shape(l,m));

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
   void get_op(const MPS_DIRECTION &dir,int site,int opnum,const QSDArray<3> &A,QSDArray<4> &E_op){

      if(dir == mpsxx::Left){

         enum {i,j,k,l,m};

         char name[100];

         sprintf(name,"scratch/Left/site_%d/%d.mpx",site,opnum);

         std::ifstream fin(name);
         boost::archive::binary_iarchive iar(fin);

         QSDArray<2> tmp2;

         iar >> tmp2;

         QSDArray<3> tmp3;

         QSDindexed_contract(1.0,tmp2,shape(l,m),A,shape(l,j,k),0.0,tmp3,shape(m,k,j));

         QSDindexed_contract(1.0,tmp3,shape(m,k,j),A.conjugate(),shape(m,i,l),0.0,E_op,shape(k,l,i,j));

      }
      else{

         enum {i,j,k,l,m};

         char name[100];

         sprintf(name,"scratch/Right/site_%d/%d.mpx",site,opnum);

         std::ifstream fin(name);
         boost::archive::binary_iarchive iar(fin);

         QSDArray<2> tmp2;

         iar >> tmp2;

         QSDArray<3> tmp3;

         QSDindexed_contract(1.0,A,shape(k,j,l),tmp2,shape(l,m),0.0,tmp3,shape(k,j,m));

         QSDindexed_contract(1.0,A.conjugate(),shape(l,i,m),tmp3,shape(k,j,m),0.0,E_op,shape(k,l,i,j));

      }

   }

   /**
    * print the contraction of two QSDArrays A on site 'site'
    */
   void print_op(const MPS_DIRECTION &dir,int site,int opnum,const QSDArray<2> &op,const QSDArray<4> &E_op){

      enum {i,j,k,l,m};

      QSDArray<2> E;
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
   void read(const MPS_DIRECTION &dir,int site,int opnum,QSDArray<2> &E){

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
   void save(const MPS_DIRECTION &dir,int site,int opnum,const QSDArray<2> &E){

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
