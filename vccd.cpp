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
      void solve(DArray<4> &t,const MPO<Q> &qc,const MPS<Q> &hf,const std::vector<double> &e,int D,double relax){

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

         int flag = 0;

         if(fabs(relax) < 1.0e-10)//line search
            flag = 1;

         double step;

         if(flag == 1)
            step = 0.2;
         else
            step = relax;

         while(convergence > 1.0e-10){

            print_dim(eTA);

            vccd::gradient(t,qc,E,eTA,eTA,list,grad,true);

            //scale the gradient: grad/M
            for(int i = 0;i < no;++i)
               for(int j = 0;j < no;++j)
                  for(int a = 0;a < nv;++a)
                     for(int b = 0;b < nv;++b)
                        grad(i,j,a,b) = grad(i,j,a,b)/M(i,j,a,b);

            if(flag == 1)
               step = line_search(qc,hf,t,grad,step,D);
            else
               step = relax;

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

         double fa = line_search_func(a,t,dir,qc,hf,D);
         double fc = line_search_func(c,t,dir,qc,hf,D);

         cout << endl;
         cout << "middle point lower than 0 point search" << endl;
         cout << endl;

         //make sure the point is lower than the 0 point
         while(fc > fa){

            c *= 0.5;

            fc = line_search_func(c,t,dir,qc,hf,D);

            cout << c << "\t" << fa << "\t" << fc << endl;

         }

         cout << endl;
         cout << "end point higher than middle point search" << endl;
         cout << endl;

         double b = c*(1.0 + phi);
         double fb = line_search_func(b,t,dir,qc,hf,D);

         //now find a b with fb which is higher than fc
         while(fc > fb){

            cout << 0 << "\t" << c << "\t" << b << "\t" << fa << "\t" << fc << "\t" << fb << endl;

            c = b;
            fc = fb;

            b = c * (1.0 + phi);

            fb = line_search_func(b,t,dir,qc,hf,D);

         }

         double d = a + b - c;
         double fd;

         cout << endl;
         cout << "actual minima search" << endl;
         cout << endl;

         while(b - a > 1.0e-4){

            fd = line_search_func(d,t,dir,qc,hf,D);

            cout << a << "\t" << c << "\t" << b << "\t" << fc << endl;

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

         return c;

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

   /**
    * optimize the cluster amplitudes t using a conjugate gradient algorithm
    * @param t input cluster amplitudes
    * @param qc the quantum-chemical hamiltonian
    * @param hf the hartree-fock reference state
    * @param e vector containing the hartree-fock energies, used to precondition the gradient
    * @param D cutoff dimension of the MPS
    */
   template<class Q>
      void conjugate_gradient(DArray<4> &t,const MPO<Q> &qc,const MPS<Q> &hf,const std::vector<double> &e,int D,int order){

         int no = t.shape(0);
         int nv = t.shape(2);

         DArray<4> M(no,no,nv,nv);

         for(int i = 0;i < no;++i)
            for(int j = 0;j < no;++j)
               for(int a = 0;a < nv;++a)
                  for(int b = 0;b < nv;++b)
                     M(i,j,a,b) = e[i] + e[j] - e[a + no] - e[b + no];

         MPO<Q> T = T2<Q>(t,false);
         compress(T,mpsxx::Right,0);
         compress(T,mpsxx::Left,0);

         MPS<Q> wccd = exp(T,hf,no,D);
         normalize(wccd);

         double E = inprod(mpsxx::Left,wccd,qc,wccd);

         cout << 0 << "\t" << E << endl;

         DArray<4> grad_1(no,no,nv,nv);
         DArray<4> grad_0(no,no,nv,nv);

         T2_2_mpo list(no,nv);

         gradient(t,qc,E,wccd,wccd,list,grad_1,true);

         //precondition the gradient
         for(int i = 0;i < no;++i)
            for(int j = 0;j < no;++j)
               for(int a = 0;a < nv;++a)
                  for(int b = 0;b < nv;++b)
                     grad_1(i,j,a,b) = grad_1(i,j,a,b)/M(i,j,a,b);

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
         double step = 0.2;

         while(convergence > 1.0e-10){

            ++iter;

            //stepsize
            step = line_search(qc,hf,t,dir,step,D);

            convergence = Ddot(grad_1,grad_1) * step;

            //update t
            Daxpy(step,dir,t);

            T = T2<Q>(t,false);
            compress(T,mpsxx::Right,0);
            compress(T,mpsxx::Left,0);

            //update the wavefunction
            wccd = exp(T,hf,no,D);
            normalize(wccd);

            cout << endl;
            print_dim(wccd);
            cout << endl;

            E = inprod(mpsxx::Left,wccd,qc,wccd);

            cout << endl;
            cout << step << "\t" << iter << "\t" << convergence << "\t" << E << endl;
            cout << endl;

            //backup the gradient
            grad_0 = grad_1;

            //calculate the gradient for the new T
            gradient(t,qc,E,wccd,wccd,list,grad_1,true);

            //precondition the gradient
            for(int i = 0;i < no;++i)
               for(int j = 0;j < no;++j)
                  for(int a = 0;a < nv;++a)
                     for(int b = 0;b < nv;++b)
                        grad_1(i,j,a,b) = grad_1(i,j,a,b)/M(i,j,a,b);

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

   template void gradient<Quantum>(const DArray<4> &,const MPO<Quantum> &,double E,const MPS<Quantum> &tccd,const MPS<Quantum> &wccd,const T2_2_mpo &list,DArray<4> &grad,bool merged);
   template void solve<Quantum>(DArray<4> &t,const MPO<Quantum> &qc,const MPS<Quantum> &hf,const std::vector<double> &,int,double);
   template double line_search<Quantum>(const MPO<Quantum> &qc,const MPS<Quantum> &hf,const DArray<4> &t,const DArray<4> &dir,double guess,int D);
   template double line_search_func<Quantum>(double a,const DArray<4> &t,const DArray<4> &dir,const MPO<Quantum> &qc,const MPS<Quantum> &hf,int D);
   template void conjugate_gradient<Quantum>(DArray<4> &t,const MPO<Quantum> &qc,const MPS<Quantum> &hf,const std::vector<double> &e,int D,int order);


}
