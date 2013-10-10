#include <iostream>
#include <iomanip>
#include <fstream>

using std::cout;
using std::endl;
using std::ostream;
using std::ofstream;

#include "include.h"

using namespace mpsxx;
using namespace btas;

namespace vccd {

   /**
    * construct the gradient of the energy for ccd
    */
   template<class Q>
      void gradient(const DArray<4> &t,const MPO<Q> &qcham,double E,const MPS<Q> &tccd,const MPS<Q> &wccd,const T2_2_mpo& list,DArray<4> &grad,bool merged){

         int no = grad.shape(0);//number of occupied orbitals
         int nv = grad.shape(2);//number of virtual orbitals

         MPO<Q> T = T2<Q>(t,merged);

         MPO<Q> rolH = ro::construct(mpsxx::Left,tccd,T,qcham,wccd);
         MPO<Q> rorH = ro::construct(mpsxx::Right,tccd,T,qcham,wccd);

         MPO<Q> mpogrH = grad::construct(rolH,rorH,tccd,qcham,wccd);

         MPS<Q> roln = ro::construct(mpsxx::Left,tccd,T,wccd);
         MPS<Q> rorn = ro::construct(mpsxx::Right,tccd,T,wccd);

         MPO<Q> mpogrn = grad::construct(roln,rorn,tccd,wccd);

         for(int i = 0;i < no;++i){

            for(int a = 0;a < nv;++a)
               for(int b = 0;b < nv;++b)
                  grad(i,i,a,b) = 2.0 * ( list.get(mpogrH,i,i,b,a,merged) + list.get(mpogrH,i,i,a,b,merged) - E * (list.get(mpogrn,i,i,b,a,merged) + list.get(mpogrn,i,i,a,b,merged)) );

            for(int j = i + 1;j < no;++j)
               for(int a = 0;a < nv;++a)
                  for(int b = 0;b < nv;++b)
                     grad(i,j,a,b) = 2.0 * ( list.get(mpogrH,i,j,a,b,merged) - E * list.get(mpogrn,i,j,a,b,merged) );

         }

         //symmetrize
         for(int i = 0;i < no;++i)
            for(int j = i + 1;j < no;++j)
               for(int a = 0;a < nv;++a)
                  for(int b = 0;b < nv;++b)
                     grad(j,i,a,b) = grad(i,j,a,b);

      }

   /**
    * solve the problem self consistently/ which in this case is more or less the same as a modified steepest descent algorithm
    */
   template<class Q>
      void solve(DArray<4> &t,const MPO<Q> &qc,const MPS<Q> &hf,const std::vector<double> &e,int D){

         int no = t.shape(0);
         int nv = t.shape(2);

         DArray<4> M(no,no,nv,nv);

         for(int i = 0;i < no;++i)
            for(int j = 0;j < no;++j)
               for(int a = 0;a < nv;++a)
                  for(int b = 0;b < nv;++b)
                     M(i,j,a,b) = e[i] + e[j] - e[a + no] - e[b + no];

         MPO<Quantum> T = T2<Quantum>(t,false);
         compress(T,mpsxx::Left,0);
         compress(T,mpsxx::Right,0);

         MPS<Quantum> eTA = exp(T,hf,no,D);
         normalize(eTA);

         double E = inprod(mpsxx::Left,eTA,qc,eTA);

         cout << 0 << "\t" << E << endl;

         DArray<4> grad(no,no,nv,nv);
         T2_2_mpo list(no,nv);

         double convergence = 1.0;

         int iter = 1;

         double step = 0.2;

         while(convergence > 1.0e-10){

            print_dim(eTA);

            vccd::gradient(t,qc,E,eTA,eTA,list,grad,true);

            //scale the gradient: grad/M
            for(int i = 0;i < no;++i)
               for(int j = 0;j < no;++j)
                  for(int a = 0;a < nv;++a)
                     for(int b = 0;b < nv;++b)
                        grad(i,j,a,b) = grad(i,j,a,b)/M(i,j,a,b);

            step = line_search(qc,hf,t,grad,step,D);

            Daxpy(step,grad,t);

            convergence = Ddot(grad,grad)*step;

            T = T2<Quantum>(t,false);
            compress(T,mpsxx::Left,0);
            compress(T,mpsxx::Right,0);

            eTA = exp(T,hf,no,D);
            normalize(eTA);

            E = inprod(mpsxx::Left,eTA,qc,eTA);

            cout << iter << "\t" << convergence << "\t" << E << endl;

            ++iter;

         }

      }

   /**
    * find the minimum of the function in the direction dir
    */
   template<class Q>
      double line_search(const MPO<Q> &qc,const MPS<Q> &hf,const DArray<4> &t,const DArray<4> &dir,double guess,int D){

         double phi = 0.5 * (1.0 + std::sqrt(5.0));

         //first bracket interval in which the minimum lies
         double a = 0.0;
         double c = guess;
         double b = c*(1.0 + phi);

         double fb = line_search_func(b,t,dir,qc,hf,D);
         double fc = line_search_func(c,t,dir,qc,hf,D);

         while(fc > fb){

            c = b;
            fc = fb;

            b = c * (1.0 + phi);

            fb = line_search_func(b,t,dir,qc,hf,D);

         }

         double d = a + b - c;
         double fd;

         while(b - a > 1.0e-4){

            fd = line_search_func(d,t,dir,qc,hf,D);

            if(c < d){

               if(fd < fc){

                  a = c;

                  c = d;
                  fc = fd;

               }
               else//fd > fc
                  b = d;

            }
            else{

               if(fd < fc){

                  b = c;

                  c = d;
                  fc = fd;

               }
               else//fd > fc
                  a = d;

            }

            //update d
            d = a + b - c;

         }

         return d;

      }

   /**
    * find the minimum of the function in the direction dir
    */
   template<class Q>
      double line_search_func(double a,const DArray<4> &t,const DArray<4> &dir,const MPO<Q> &qc,const MPS<Q> &hf,int D){

         DArray<4> newt(t);
         Daxpy(a,dir,newt);

         int no = newt.shape(0);

         int L = qc.size();

         MPO<Q> T2op = T2<Quantum>(newt,false);
         compress(T2op,mpsxx::Right,0);
         compress(T2op,mpsxx::Left,0);

         MPS<Q> eTA = exp(T2op,hf,no,D);
         normalize(eTA);

         return inprod(mpsxx::Left,eTA,qc,eTA);

      }

   template void gradient<Quantum>(const DArray<4> &,const MPO<Quantum> &,double E,const MPS<Quantum> &tccd,const MPS<Quantum> &wccd,const T2_2_mpo &list,DArray<4> &grad,bool merged);
   template void solve<Quantum>(DArray<4> &t,const MPO<Quantum> &qc,const MPS<Quantum> &hf,const std::vector<double> &,int);
   template double line_search<Quantum>(const MPO<Quantum> &qc,const MPS<Quantum> &hf,const DArray<4> &t,const DArray<4> &dir,double guess,int D);
   template double line_search_func<Quantum>(double a,const DArray<4> &t,const DArray<4> &dir,const MPO<Quantum> &qc,const MPS<Quantum> &hf,int D);

}
