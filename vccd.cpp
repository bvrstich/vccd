#include <iostream>
#include <iomanip>

using std::cout;
using std::endl;
using std::ostream;

#include "include.h"

using namespace mpsxx;
using namespace btas;

namespace vccd {

   /**
    * construct the gradient of the energy for ccd
    */
   template<class Q>
      void gradient(const MPO<Q> &qcham,const MPS<Q> &wccd,DArray<4> &grad){
/*
         double norm_ccd = dot(mpsxx::Left,wccd,wccd);
         MPS<Q> Hccd = qcham*wccd;
         compress(Hccd,mpsxx::Right,0);
         compress(Hccd,mpsxx::Left,0);

         double E_ccd = dot(mpsxx::Left,wccd,Hccd)/norm_ccd;

         int no = grad.shape(0);//number of occupied orbitals
         int nv = grad.shape(2);//number of virtual orbitals

         int L = no + nv;

         MPO<Q> Eabij;
         MPO<Q> prod;
         MPS<Q> tmp;

         for(int i = 0;i < no;++i)
            for(int j = 0;j < no;++j)
               for(int a = 0;a < nv;++a)
                  for(int b = 0;b < nv;++b){

                     grad(i,j,a,b) = 0.0;

                     Eabij = E<Quantum>(L,a+no,b+no,i,j,1.0);
                     tmp = Eabij*wccd;

                     grad(i,j,a,b) = 2.0 * (dot(mpsxx::Left,tmp,Hccd) - E_ccd * dot(mpsxx::Left,tmp,wccd))/norm_ccd;

                  }
*/

         int no = grad.shape(0);//number of occupied orbitals
         int nv = grad.shape(2);//number of virtual orbitals

         int L = no + nv;

         MPO<Q> Eabij;
         MPO<Q> prod;
         MPS<Q> tmp;

         for(int i = 0;i < no;++i)
            for(int j = 0;j < no;++j)
               for(int a = 0;a < nv;++a)
                  for(int b = 0;b < nv;++b){

                     grad(i,j,a,b) = 0.0;

                     Eabij = E<Quantum>(L,a+no,b+no,i,j,1.0);
                     tmp = Eabij*wccd;

                     grad(i,j,a,b) = dot(mpsxx::Left,tmp,wccd);

                     cout << i << "\t" << j << "\t" << a << "\t" << b << "\t" << grad(i,j,a,b) << endl;

                  }

      }

   /**
    * construct the gradient of the energy for ccd
    */
   template<class Q>
      void gradient_new(const DArray<4> &t,const MPO<Q> &qcham,const MPS<Q> &wccd,const T2_2_mpo& list,DArray<4> &grad){

         int no = grad.shape(0);//number of occupied orbitals
         int nv = grad.shape(2);//number of virtual orbitals

         MPO<Q> T = T2<Quantum>(t,false);

         MPS<Quantum> rol = ro::construct(mpsxx::Left,wccd,T,wccd);
         MPS<Quantum> ror = ro::construct(mpsxx::Right,wccd,T,wccd);

         MPO<Quantum> mpogr = grad::construct(rol,ror,wccd,wccd);

         for(int i = 0;i < no;++i){

            for(int a = 0;a < nv;++a)
               for(int b = 0;b < nv;++b)
                  grad(i,i,a,b) = list.get(mpogr,i,i,b,a) + list.get(mpogr,i,i,a,b);

            for(int j = i + 1;j < no;++j)
               for(int a = 0;a < nv;++a)
                  for(int b = 0;b < nv;++b)
                     grad(i,j,a,b) = list.get(mpogr,i,j,a,b);

         }

         //symmetrize
         for(int i = 0;i < no;++i)
            for(int j = i + 1;j < no;++j)
               for(int a = 0;a < nv;++a)
                  for(int b = 0;b < nv;++b)
                     grad(j,i,a,b) = grad(i,j,a,b);

      }


   /**
    * find the minimum of the function in the direction dir
    */
   template<class Q>
      double line_search(const MPO<Q> &qc,const MPS<Q> &hf,const DArray<4> &t,const DArray<4> &dir,double guess,const std::vector<int> &cutoff){

         double phi = 0.5 * (1.0 + std::sqrt(5.0));

         //first bracket interval in which the minimum lies
         double a = 0.0;
         double c = guess;
         double b = c*(1.0 + phi);

         double fb = line_search_func(b,t,dir,qc,hf,cutoff);
         double fc = line_search_func(c,t,dir,qc,hf,cutoff);

         while(fc > fb){

            c = b;
            fc = fb;

            b = c * (1.0 + phi);

            fb = line_search_func(b,t,dir,qc,hf,cutoff);

         }

         double d = a + b - c;
         double fd;

         while(b - a > 1.0e-4){

            fd = line_search_func(d,t,dir,qc,hf,cutoff);

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
      double line_search_func(double a,const DArray<4> &t,const DArray<4> &dir,const MPO<Q> &qc,const MPS<Q> &hf,const std::vector<int> &cutoff){

         DArray<4> newt(t);
         Daxpy(-a,dir,newt);

         MPO<Q> T2op = T2<Quantum>(newt,true);
         compress(T2op,mpsxx::Right,cutoff[0]);
         compress(T2op,mpsxx::Left,cutoff[0]);

         MPS<Q> eTA = exp(T2op,hf,cutoff);
         normalize(eTA);

         return inprod(mpsxx::Left,eTA,qc,eTA);

      }

   template<class Q>
      void steepest_descent(DArray<4> &t,const MPO<Q> &qc,const MPS<Q> &hf,const std::vector<int> &cutoff){

         int no = t.shape(0);
         int nv = t.shape(2);

         MPO<Quantum> T2op = T2<Quantum>(t,true);
         compress(T2op,mpsxx::Right,0);
         compress(T2op,mpsxx::Left,0);

         MPS<Quantum> eTA = exp(T2op,hf,cutoff);
         normalize(eTA);

         DArray<4> grad(no,no,nv,nv);
         gradient(qc,eTA,grad);

         T2_2_mpo list(no,nv);

         gradient_new(t,qc,eTA,list,grad);

         /*
            double step = line_search(qc,hf,t,grad,0.4,cutoff);

            double convergence = 1.0;

            while(convergence > 1.0e-5){

            Daxpy(-step,grad,t);

            T2op = T2<Quantum>(t,true);
            compress(T2op,mpsxx::Right,0);
            compress(T2op,mpsxx::Left,0);

            eTA = exp(T2op,hf,cutoff);
            normalize(eTA);

            convergence = Ddot(grad,grad);

            cout << convergence << "\t" << inprod(mpsxx::Left,eTA,qc,eTA) << endl;

            gradient(qc,eTA,grad);
            step = line_search(qc,hf,t,grad,0.2,cutoff);

            }
          */
      }

   template void gradient<Quantum>(const MPO<Quantum> &,const MPS<Quantum> &wccd,DArray<4> &grad);
   template void gradient_new<Quantum>(const DArray<4> &,const MPO<Quantum> &,const MPS<Quantum> &wccd,const T2_2_mpo &list,DArray<4> &grad);
   template double line_search<Quantum>(const MPO<Quantum> &,const MPS<Quantum> &,const DArray<4> &,const DArray<4> &,double,const std::vector<int> &cutoff);
   template double line_search_func<Quantum>(double a,const DArray<4> &t,const DArray<4> &dir,const MPO<Quantum> &qc,const MPS<Quantum> &hf,const std::vector<int> &cutoff);
   template void steepest_descent<Quantum>(DArray<4>  &,const MPO<Quantum> &qc,const MPS<Quantum> &hf,const std::vector<int> &cutoff);

}
