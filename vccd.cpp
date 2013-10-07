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
                  grad(i,i,a,b) = -2.0 * ( list.get(mpogrH,i,i,b,a,merged) + list.get(mpogrH,i,i,a,b,merged) - E * (list.get(mpogrn,i,i,b,a,merged) + list.get(mpogrn,i,i,a,b,merged)) );

            for(int j = i + 1;j < no;++j)
               for(int a = 0;a < nv;++a)
                  for(int b = 0;b < nv;++b)
                     grad(i,j,a,b) = -2.0 * ( list.get(mpogrH,i,j,a,b,merged) - E * list.get(mpogrn,i,j,a,b,merged) );

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

         while(b - a > 1.0e-10){

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

         MPO<Q> T = T2<Q>(t,false);
         compress(T,mpsxx::Right,0);
         compress(T,mpsxx::Left,0);

         eMPS emps(T,hf,cutoff);

         MPS<Q> wccd = emps.expand(hf,cutoff.size(),0);
         MPS<Q> tccd = emps.expand(hf,cutoff.size() - 1,0);

         cout << endl;
         emps.print_tot_dim(0);
         cout << endl;

         //get the norm right
         double nrm = mpsxx::nrm2(wccd);

         mpsxx::scal(1.0/nrm,wccd);
         mpsxx::scal(1.0/nrm,tccd);

         //double E = emps.eval(qc,hf);
         //double N = emps.norm();
         double E = inprod(mpsxx::Left,wccd,qc,wccd);
         double N = dot(mpsxx::Left,wccd,wccd);
         cout << E/N << endl;

         DArray<4> grad(no,no,nv,nv);

         T2_2_mpo list(no,nv);

         //gradient(t,qc,E/N,tccd,wccd,list,grad,false);
         //old_gradient(qc,wccd,grad);
         num_gradient(E/N,qc,hf,t,grad,cutoff);

         //double step = line_search(E,N,qc,hf,emps,T,0.4);
         double step  = old_line_search(qc,hf,t,grad,0.01,cutoff);

         double convergence = 1.0;

         int iter = 0;

         while(convergence > 1.0e-10){

            ++iter;

            Daxpy(step,grad,t);

            T = T2<Q>(t,false);
            compress(T,mpsxx::Right,0);
            compress(T,mpsxx::Left,0);

            //set new T
            emps.update(T,hf);

            //get the wavefunction out
            wccd = emps.expand(hf,cutoff.size(),0);
            tccd = emps.expand(hf,cutoff.size() - 1,0);

            cout << endl;
            emps.print_tot_dim(0);
            cout << endl;

            //get the norm right
            nrm = mpsxx::nrm2(wccd);

            mpsxx::scal(1.0/nrm,wccd);
            mpsxx::scal(1.0/nrm,tccd);

            normalize(wccd);

            convergence = Ddot(grad,grad) * step * step;

            //E = emps.eval(qc,hf);
            //N = emps.norm();

            E = mpsxx::inprod(mpsxx::Left,wccd,qc,wccd);
            N = mpsxx::dot(mpsxx::Left,wccd,wccd);

            cout << iter << "\t" << convergence << "\t" << E/N << endl;

            //gradient(t,qc,E/N,tccd,wccd,list,grad,false);
            //old_gradient(qc,wccd,grad);
            num_gradient(E/N,qc,hf,t,grad,cutoff);

            T = T2<Q>(grad,false);
            compress(T,mpsxx::Right,0);
            compress(T,mpsxx::Left,0);

            //step = line_search(E,N,qc,hf,emps,T,step);
            step  = old_line_search(qc,hf,t,grad,step,cutoff);

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

         cout << endl;
         emps.print_tot_dim(0);
         cout << endl;

         MPS<Q> wccd = emps.expand(hf,cutoff.size(),0);
         MPS<Q> tccd = emps.expand(hf,cutoff.size() - 1,0);

         //get the norm right
         double nrm = mpsxx::nrm2(wccd);

         mpsxx::scal(1.0/nrm,wccd);
         mpsxx::scal(1.0/nrm,tccd);

         DArray<4> grad_1(no,no,nv,nv);
         DArray<4> grad_0(no,no,nv,nv);

         T2_2_mpo list(no,nv);

         //gradient(t,qc,E/N,tccd,wccd,list,grad_1,true);
         old_gradient(qc,wccd,grad_1,cutoff);

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

         while(convergence > 1.0e-10){

            ++iter;

            //stepsize
            //step = line_search(E,N,qc,hf,emps,T,step);
            step = old_line_search(qc,hf,t,dir,step,cutoff);

            convergence = Ddot(grad_1,grad_1) * step * step;

            //update t
            Daxpy(step,dir,t);

            T = T2<Q>(t,false);
            compress(T,mpsxx::Right,0);
            compress(T,mpsxx::Left,0);

            //set new T
            emps.update(T,hf);

            cout << endl;
            emps.print_tot_dim(0);
            cout << endl;

            E = emps.eval(qc,hf);
            N = emps.norm();
            cout << iter << "\t" << convergence << "\t" << E/N << endl;

            //get the wavefunction out
            wccd = emps.expand(hf,cutoff.size(),0);
            tccd = emps.expand(hf,cutoff.size() - 1,0);

            //get the norm right
            double nrm = mpsxx::nrm2(wccd);

            mpsxx::scal(1.0/nrm,wccd);
            mpsxx::scal(1.0/nrm,tccd);

            //backup the gradient
            grad_0 = grad_1;

            //calculate the gradient for the new T
            //gradient(t,qc,E/N,tccd,wccd,list,grad_1,true);
            old_gradient(qc,wccd,grad_1,cutoff);

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

   /**
    * optimize the cluster amplitudes t using a conjugate gradient algorithm
    * @param t input cluster amplitudes
    * @param qc the quantum-chemical hamiltonian
    * @param hf the hartree-fock reference state
    * @param cutoff vector containing the number of states to be withheld for every order in the expansion
    */
   template<class Q>
      void conjugate_gradient_tmp(DArray<4> &t,const MPO<Q> &qc,const MPS<Q> &hf,const std::vector<int> &cutoff){

         int no = t.shape(0);
         int nv = t.shape(2);

         MPO<Q> T = T2<Q>(t,false);
         compress(T,mpsxx::Right,0);
         compress(T,mpsxx::Left,0);

         eMPS emps(T,hf,cutoff);

         double E = emps.eval(qc,hf);
         double N = emps.norm();
         cout << E/N << endl;

         cout << endl;
         emps.print_tot_dim(0);
         cout << endl;

         MPS<Q> wccd = emps.expand(hf,cutoff.size(),0);
         MPS<Q> tccd = emps.expand(hf,cutoff.size() - 1,0);

         //get the norm right
         double nrm = mpsxx::nrm2(wccd);

         mpsxx::scal(1.0/nrm,wccd);
         mpsxx::scal(1.0/nrm,tccd);

         DArray<4> grad_1(no,no,nv,nv);
         DArray<4> grad_0(no,no,nv,nv);

         T2_2_mpo list(no,nv);

         //gradient(t,qc,E/N,tccd,wccd,list,grad_1,true);
         num_gradient(E/N,qc,hf,t,grad_1,cutoff);

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
         double step = 0.1;

         while(convergence > 1.0e-10){

            ++iter;

            //stepsize
            //step = line_search(E,N,qc,hf,emps,T,step);
            step  = old_line_search(qc,hf,t,dir,step,cutoff);

            convergence = Ddot(grad_1,grad_1) * step * step;

            //update t
            Daxpy(step,dir,t);

            T = T2<Q>(t,false);
            compress(T,mpsxx::Right,0);
            compress(T,mpsxx::Left,0);

            //set new T
            emps.update(T,hf);

            cout << endl;
            emps.print_tot_dim(0);
            cout << endl;

            E = emps.eval(qc,hf);
            N = emps.norm();
            cout << iter << "\t" << convergence << "\t" << E/N << endl;

            //get the wavefunction out
            wccd = emps.expand(hf,cutoff.size(),0);
            tccd = emps.expand(hf,cutoff.size() - 1,0);

            //get the norm right
            double nrm = mpsxx::nrm2(wccd);

            mpsxx::scal(1.0/nrm,wccd);
            mpsxx::scal(1.0/nrm,tccd);

            //backup the gradient
            grad_0 = grad_1;

            //calculate the gradient for the new T
            //   gradient(t,qc,E/N,tccd,wccd,list,grad_1,true);
            num_gradient(E/N,qc,hf,t,grad_1,cutoff);

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


   /**
    * find the minimum of the function in the direction dir
    */
   template<class Q>
      double old_line_search(const MPO<Q> &qc,const MPS<Q> &hf,const DArray<4> &t,const DArray<4> &dir,double guess,const std::vector<int> &cutoff){

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
         Daxpy(a,dir,newt);

         int L = qc.size();

         MPO<Q> T2op = T2<Quantum>(newt,false);
         compress(T2op,mpsxx::Right,0);
         compress(T2op,mpsxx::Left,0);

         eMPS emps(T2op,hf,cutoff);
         MPS<Q> eTA = emps.expand(hf,cutoff.size(),0);
         normalize(eTA);

         return inprod(mpsxx::Left,eTA,qc,eTA);

      }

   /**
    * construct the gradient of the energy for ccd
    */
   template<class Q>
      void old_gradient(const MPO<Q> &qcham,const MPS<Q> &wccd,DArray<4> &grad,const std::vector<int> &cutoff){

         double norm_ccd = dot(mpsxx::Left,wccd,wccd);
         MPS<Q> Hccd = qcham*wccd;

         compress(Hccd,mpsxx::Left,0);
         compress(Hccd,mpsxx::Right,0);

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

                     compress(tmp,mpsxx::Left,0);
                     compress(tmp,mpsxx::Right,cutoff[0] + 1);

                     grad(i,j,a,b) = -2.0 * (dot(mpsxx::Left,tmp,Hccd) - E_ccd * dot(mpsxx::Left,tmp,wccd))/norm_ccd;

                     cout << i << "\t" << j << "\t" << a << "\t" << b << "\t"  << grad(i,j,a,b) << endl;

                  }


      }

   /**
    * construct the gradient of the energy for ccd
    */
   template<class Q>
      void num_gradient(double E,const MPO<Q> &qcham,const MPS<Q> &hf,const DArray<4> &t,DArray<4> &grad,const std::vector<int> &cutoff){

         MPO<Q> T;

         int no = grad.shape(0);//number of occupied orbitals
         int nv = grad.shape(2);//number of virtual orbitals

         int L = no + nv;

         DArray<4> tmp;
         MPS<Quantum> eTA;

         double h = 0.0001;

         for(int i = 0;i < no;++i)
            for(int j = i;j < no;++j)
               for(int a = 0;a < nv;++a)
                  for(int b = 0;b < nv;++b){

                     grad(i,j,a,b) = 0.0;

                     tmp = t;
                     tmp(i,j,a,b) += h;

                     if(i != j && b != a)
                        tmp(j,i,b,a) += h;

                     T = T2<Quantum>(tmp,false);
                     compress(T,mpsxx::Left,0);
                     compress(T,mpsxx::Right,0);

                     eTA = exp(T,hf,cutoff);

                     double nrm = mpsxx::dot(mpsxx::Left,eTA,eTA);

                     grad(i,j,a,b) = (E - inprod(mpsxx::Left,eTA,qcham,eTA)/nrm)/h;

                     cout << i << "\t" << j << "\t" << a << "\t" << b << "\t" << grad(i,j,a,b) << endl;

                  }

         for(int i = 0;i < no;++i)
            for(int j = i + 1;j < no;++j)
               for(int a = 0;a < nv;++a)
                  for(int b = 0;b < nv;++b)
                     grad(j,i,a,b) = grad(i,j,b,a);


      }

   template void gradient<Quantum>(const DArray<4> &,const MPO<Quantum> &,double E,const MPS<Quantum> &tccd,const MPS<Quantum> &wccd,const T2_2_mpo &list,DArray<4> &grad,bool merged);
   template double line_search<Quantum>(double,double,const MPO<Quantum> &,const MPS<Quantum> &,const eMPS &,const MPO<Quantum> &,double );
   template void steepest_descent<Quantum>(DArray<4>  &,const MPO<Quantum> &qc,const MPS<Quantum> &hf,const std::vector<int> &cutoff);
   template void conjugate_gradient<Quantum>(DArray<4>  &,const MPO<Quantum> &qc,const MPS<Quantum> &hf,const std::vector<int> &cutoff);
   template void conjugate_gradient_tmp<Quantum>(DArray<4>  &,const MPO<Quantum> &qc,const MPS<Quantum> &hf,const std::vector<int> &cutoff);
   template double old_line_search(const MPO<Quantum> &qc,const MPS<Quantum> &hf,const DArray<4> &t,const DArray<4> &dir,double guess,const std::vector<int> &cutoff);
   template double line_search_func(double a,const DArray<4> &t,const DArray<4> &dir,const MPO<Quantum> &qc,const MPS<Quantum> &hf,const std::vector<int> &cutoff);
   template void old_gradient(const MPO<Quantum> &qcham,const MPS<Quantum> &wccd,DArray<4> &grad,const std::vector<int> &);
   template void num_gradient(double,const MPO<Quantum> &qcham,const MPS<Quantum> &hf,const DArray<4> &t,DArray<4> &grad,const std::vector<int> &);

}
