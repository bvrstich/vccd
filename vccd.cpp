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
      void gradient(const DArray<4> &t,const MPO<Q> &qcham,double E,const MPS<Q> &wccd,const T2_2_mpo& list,DArray<4> &grad){

         int no = grad.shape(0);//number of occupied orbitals
         int nv = grad.shape(2);//number of virtual orbitals

         MPS<Q> Hccd = qcham*wccd;
         compress(Hccd,mpsxx::Right,0);
         compress(Hccd,mpsxx::Left,0);

         MPO<Q> T = T2<Quantum>(t,false);

         MPS<Quantum> rolH = ro::construct(mpsxx::Left,wccd,T,Hccd);
         MPS<Quantum> rorH = ro::construct(mpsxx::Right,wccd,T,Hccd);

         MPO<Quantum> mpogrH = grad::construct(rolH,rorH,wccd,Hccd);

         MPS<Quantum> roln = ro::construct(mpsxx::Left,wccd,T,wccd);
         MPS<Quantum> rorn = ro::construct(mpsxx::Right,wccd,T,wccd);

         MPO<Quantum> mpogrn = grad::construct(roln,rorn,wccd,wccd);

         for(int i = 0;i < no;++i){

            for(int a = 0;a < nv;++a)
               for(int b = 0;b < nv;++b)
                  grad(i,i,a,b) = 2.0 * ( list.get(mpogrH,i,i,b,a) + list.get(mpogrH,i,i,a,b) - E * (list.get(mpogrn,i,i,b,a) + list.get(mpogrn,i,i,a,b)) );

            for(int j = i + 1;j < no;++j)
               for(int a = 0;a < nv;++a)
                  for(int b = 0;b < nv;++b)
                     grad(i,j,a,b) = 2.0 * ( list.get(mpogrH,i,j,a,b) - E * list.get(mpogrn,i,j,a,b) );

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
      double new_line_search(double E,double N,const MPO<Q> &qc,const MPS<Q> &hf,const eMPS &ccd,const MPO<Q> &dir){

         e_eMPS e_emps(dir,ccd,hf);

         int matdim = ccd.gcutoff().size() + 1;

         DArray<2> Emat(matdim,matdim);
         DArray<2> Nmat(matdim,matdim);

         Emat(0,0) = E;
         Nmat(0,0) = N;

         ls::construct(Emat,Nmat,qc,hf,ccd,e_emps);

         cout << Emat << endl;
         cout << Nmat << endl;
/*
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
*/
         return 0.0;
      }



   /**
    * find the minimum of the function in the direction dir
    */
   template<class Q>
      double line_search_func(double a,const DArray<4> &t,const DArray<4> &dir,const MPO<Q> &qc,const MPS<Q> &hf,const std::vector<int> &cutoff){

         DArray<4> newt(t);
         Daxpy(-a,dir,newt);

         MPO<Q> T = T2<Quantum>(newt,true);
         compress(T,mpsxx::Right,cutoff[0]);
         compress(T,mpsxx::Left,cutoff[0]);

         MPS<Q> eTA = exp(T,hf,cutoff);
         normalize(eTA);

         return inprod(mpsxx::Left,eTA,qc,eTA);

      }

   template<class Q>
      void steepest_descent(DArray<4> &t,const MPO<Q> &qc,const MPS<Q> &hf,const std::vector<int> &cutoff){

         int no = t.shape(0);
         int nv = t.shape(2);

         MPO<Quantum> T = T2<Quantum>(t,true);
         compress(T,mpsxx::Right,0);
         compress(T,mpsxx::Left,0);

         eMPS emps(T,hf,cutoff);

         MPS<Quantum> wccd = emps.expand(hf,cutoff[0]);
         normalize(wccd);

         double E = emps.eval(qc,hf);
         double N = emps.norm();
         cout << E/N << endl;

         DArray<4> grad(no,no,nv,nv);

         T2_2_mpo list(no,nv);

         gradient(t,qc,E/N,wccd,list,grad);

         T = T2<Quantum>(grad,false);
         compress(T,mpsxx::Right,0);
         compress(T,mpsxx::Left,0);

         double step = line_search(qc,hf,t,grad,0.4,cutoff);
         cout << step << endl;

         step = new_line_search(E,N,qc,hf,emps,T);
         cout << step << endl;
         /*
            double convergence = 1.0;

         //         while(convergence > 1.0e-5){

         Daxpy(-step,grad,t);

         T = T2<Quantum>(t,true);
            compress(T,mpsxx::Right,0);
            compress(T,mpsxx::Left,0);

            //set new T
            emps.update(T);

            //get the wavefunction out
            wccd = emps.expand(cutoff[0]);
            normalize(wccd);

            convergence = Ddot(grad,grad);
            E = emps.eval(qc)/emps.norm();

            cout << convergence << "\t" << E << endl;

            gradient(t,qc,E,wccd,list,grad);

            step = line_search(qc,hf,t,grad,0.2,cutoff);
            cout << step << endl;

            step = new_line_search(qc,emps,grad);
            cout << step << endl;

 //        }
*/
      }

   template void gradient<Quantum>(const DArray<4> &,const MPO<Quantum> &,double E,const MPS<Quantum> &wccd,const T2_2_mpo &list,DArray<4> &grad);
   template double line_search<Quantum>(const MPO<Quantum> &,const MPS<Quantum> &,const DArray<4> &,const DArray<4> &,double,const std::vector<int> &cutoff);
   template double line_search_func<Quantum>(double a,const DArray<4> &t,const DArray<4> &dir,const MPO<Quantum> &qc,const MPS<Quantum> &hf,const std::vector<int> &cutoff);
   template double new_line_search<Quantum>(double,double,const MPO<Quantum> &,const MPS<Quantum> &,const eMPS &,const MPO<Quantum> &);
   template void steepest_descent<Quantum>(DArray<4>  &,const MPO<Quantum> &qc,const MPS<Quantum> &hf,const std::vector<int> &cutoff);

}
