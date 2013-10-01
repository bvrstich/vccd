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
    * construct minus the gradient of the energy for ccd
    */
   template<class Q>
      void gradient(const DArray<4> &t,const MPO<Q> &qcham,double E,const MPS<Q> &tccd,const MPS<Q> &wccd,const T2_2_mpo& list,DArray<4> &grad){

         int no = grad.shape(0);//number of occupied orbitals
         int nv = grad.shape(2);//number of virtual orbitals

         MPO<Q> T = T2<Q>(t,false);

         MPO<Q> rolH = ro::construct(mpsxx::Left,tccd,T,qcham,wccd);
         MPO<Q> rorH = ro::construct(mpsxx::Right,tccd,T,qcham,wccd);

         MPO<Q> mpogrH = grad::construct(rolH,rorH,tccd,qcham,wccd);

         MPS<Q> roln = ro::construct(mpsxx::Left,tccd,T,wccd);
         MPS<Q> rorn = ro::construct(mpsxx::Right,tccd,T,wccd);

         MPO<Q> mpogrn = grad::construct(roln,rorn,tccd,wccd);

         for(int i = 0;i < no;++i){

            for(int a = 0;a < nv;++a)
               for(int b = 0;b < nv;++b)
                  grad(i,i,a,b) = -2.0 * ( list.get(mpogrH,i,i,b,a) + list.get(mpogrH,i,i,a,b) - E * (list.get(mpogrn,i,i,b,a) + list.get(mpogrn,i,i,a,b)) );

            for(int j = i + 1;j < no;++j)
               for(int a = 0;a < nv;++a)
                  for(int b = 0;b < nv;++b)
                     grad(i,j,a,b) = -2.0 * ( list.get(mpogrH,i,j,a,b) - E * list.get(mpogrn,i,j,a,b) );

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
      double line_search(double E,double N,const MPO<Q> &qc,const MPS<Q> &hf,const eMPS &ccd,const MPO<Q> &dir,double guess){

         e_eMPS e_emps(dir,ccd,hf);

         int matdim = ccd.gcutoff().size() + 1;

         DArray<2> Emat(matdim,matdim);
         DArray<2> Nmat(matdim,matdim);

         Emat(0,0) = E;
         Nmat(0,0) = N;

         ls::construct(Emat,Nmat,qc,hf,ccd,e_emps);

         double phi = 0.5 * (1.0 + std::sqrt(5.0));

         //first bracket interval in which the minimum lies
         double a = 0.0;
         double c = guess;
         double b = c*(1.0 + phi);

         double fb = ls::eval(Emat,Nmat,b);
         double fc = ls::eval(Emat,Nmat,c);

         while(fc > fb){

            c = b;
            fc = fb;

            b = c * (1.0 + phi);

            fb = ls::eval(Emat,Nmat,b);

         }

         double d = a + b - c;
         double fd;

         while(b - a > 1.0e-7){

            fd = ls::eval(Emat,Nmat,d);

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
    * optimize the cluster amplitudes t using a simplest steepest descent algorithm
    * @param t input cluster amplitudes
    * @param qc the quantum-chemical hamiltonian
    * @param hf the hartree-fock reference state
    * @param cutoff vector containing the number of states to be withheld for every order in the expansion
    */
   template<class Q>
      void steepest_descent(DArray<4> &t,const MPO<Q> &qc,const MPS<Q> &hf,const std::vector<int> &cutoff){

         int no = t.shape(0);
         int nv = t.shape(2);

         MPO<Q> T = T2<Q>(t,true);
         compress(T,mpsxx::Right,0);
         compress(T,mpsxx::Left,0);

         eMPS emps(T,hf,cutoff);

         MPS<Q> wccd = emps.expand(hf,cutoff.size(),cutoff[0]);
         MPS<Q> tccd = emps.expand(hf,1,cutoff[0]);

         //get the norm right
         double nrm = mpsxx::nrm2(wccd);

         mpsxx::scal(1.0/nrm,wccd);
         mpsxx::scal(1.0/nrm,tccd);

         double E = emps.eval(qc,hf);
         double N = emps.norm();
         cout << E/N << endl;

         DArray<4> grad(no,no,nv,nv);

         T2_2_mpo list(no,nv);

         gradient(t,qc,E/N,tccd,wccd,list,grad);

         T = T2<Q>(grad,false);
         compress(T,mpsxx::Right,0);
         compress(T,mpsxx::Left,0);

         double step = line_search(E,N,qc,hf,emps,T,0.4);

         double convergence = 1.0;

         while(convergence > 1.0e-5){

            Daxpy(step,grad,t);

            T = T2<Q>(t,false);
            compress(T,mpsxx::Right,0);
            compress(T,mpsxx::Left,0);

            //set new T
            emps.update(T,hf);

            //get the wavefunction out
            wccd = emps.expand(hf,cutoff.size(),cutoff[0]);
            tccd = emps.expand(hf,1,cutoff[0]);

            //get the norm right
            double nrm = mpsxx::nrm2(wccd);

            mpsxx::scal(1.0/nrm,wccd);
            mpsxx::scal(1.0/nrm,tccd);

            normalize(wccd);

            convergence = Ddot(grad,grad);
            E = emps.eval(qc,hf);
            N = emps.norm();

            cout << convergence << "\t" << E/N << endl;

            gradient(t,qc,E/N,tccd,wccd,list,grad);

            T = T2<Q>(grad,false);
            compress(T,mpsxx::Right,0);
            compress(T,mpsxx::Left,0);

            step = line_search(E,N,qc,hf,emps,T,step);

         }

      }

   /**
    * optimize the cluster amplitudes t using a conjugate gradient algorithm
    * @param t input cluster amplitudes
    * @param qc the quantum-chemical hamiltonian
    * @param hf the hartree-fock reference state
    * @param cutoff vector containing the number of states to be withheld for every order in the expansion
    */
   template<class Q>
      void conjugate_gradient(DArray<4> &t,const MPO<Q> &qc,const MPS<Q> &hf,const std::vector<int> &cutoff){

         int no = t.shape(0);
         int nv = t.shape(2);

         MPO<Q> T = T2<Q>(t,true);
         compress(T,mpsxx::Right,0);
         compress(T,mpsxx::Left,0);

         eMPS emps(T,hf,cutoff);

         double E = emps.eval(qc,hf);
         double N = emps.norm();
         cout << E/N << endl;

         MPS<Q> wccd = emps.expand(hf,cutoff.size(),cutoff[0]);
         MPS<Q> tccd = emps.expand(hf,1,cutoff[0]);

         //get the norm right
         double nrm = mpsxx::nrm2(wccd);

         mpsxx::scal(1.0/nrm,wccd);
         mpsxx::scal(1.0/nrm,tccd);

         DArray<4> grad_1(no,no,nv,nv);
         DArray<4> grad_0(no,no,nv,nv);

         T2_2_mpo list(no,nv);

         gradient(t,qc,E/N,tccd,wccd,list,grad_1);

         T = T2<Q>(grad_1,false);
         compress(T,mpsxx::Right,0);
         compress(T,mpsxx::Left,0);

         //conjugated gradient search direction
         DArray<4> dir(grad_1);

         double rnorm_0,rnorm_1;

         //norm of the gradient
         rnorm_0 = Ddot(grad_1,grad_1);

         int iter = 0;

         double convergence = 1.0;
         double step = 0.4;

         while(convergence > 1.0e-7){

            ++iter;

            //stepsize
            step = line_search(E,N,qc,hf,emps,T,step);

            convergence = Ddot(grad_1,grad_1) * step * step;

            //update t
            Daxpy(step,dir,t);

            T = T2<Q>(t,false);
            compress(T,mpsxx::Right,0);
            compress(T,mpsxx::Left,0);

            //set new T
            emps.update(T,hf);

            E = emps.eval(qc,hf);
            N = emps.norm();
            cout << iter << "\t" << convergence << "\t" << E/N << endl;

            //get the wavefunction out
            wccd = emps.expand(hf,cutoff.size(),cutoff[0]);
            tccd = emps.expand(hf,1,cutoff[0]);

            //get the norm right
            double nrm = mpsxx::nrm2(wccd);

            mpsxx::scal(1.0/nrm,wccd);
            mpsxx::scal(1.0/nrm,tccd);

            //backup the gradient
            grad_0 = grad_1;

            //calculate the gradient for the new T
            gradient(t,qc,E/N,tccd,wccd,list,grad_1);

            rnorm_1 = Ddot(grad_1,grad_1);

            //polak-ribiery update
            double beta = (rnorm_1 - Ddot(grad_0,grad_1))/rnorm_0;

            if(beta < 0.0)
               beta = 0.0;

            //update the direction
            Dscal(beta,dir);
            Daxpy(1.0,grad_1,dir);

            rnorm_0 = rnorm_1;

            //construct the T for the line search
            T = T2<Q>(dir,false);
            compress(T,mpsxx::Right,0);
            compress(T,mpsxx::Left,0);

         }

      }

   template void gradient<Quantum>(const DArray<4> &,const MPO<Quantum> &,double E,const MPS<Quantum> &tccd,const MPS<Quantum> &wccd,const T2_2_mpo &list,DArray<4> &grad);
   template double line_search<Quantum>(double,double,const MPO<Quantum> &,const MPS<Quantum> &,const eMPS &,const MPO<Quantum> &,double );
   template void steepest_descent<Quantum>(DArray<4>  &,const MPO<Quantum> &qc,const MPS<Quantum> &hf,const std::vector<int> &cutoff);
   template void conjugate_gradient<Quantum>(DArray<4>  &,const MPO<Quantum> &qc,const MPS<Quantum> &hf,const std::vector<int> &cutoff);

}
