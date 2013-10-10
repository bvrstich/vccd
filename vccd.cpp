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

        MPS<Quantum> eTA = exp(T,hf,2,D);
        normalize(eTA);

        double E = inprod(mpsxx::Left,eTA,qc,eTA);

        cout << 0 << "\t" << E << endl;

        DArray<4> grad(no,no,nv,nv);
        T2_2_mpo list(no,nv);

        double convergence = 1.0;

        int iter = 1;

        while(convergence > 1.0e-10){

           print_dim(eTA);

           vccd::gradient(t,qc,E,eTA,eTA,list,grad,true);

           //new t is old + grad/M
           for(int i = 0;i < no;++i)
              for(int j = 0;j < no;++j)
                 for(int a = 0;a < nv;++a)
                    for(int b = 0;b < nv;++b)
                       t(i,j,a,b) += relax * grad(i,j,a,b)/M(i,j,a,b);

           convergence = Ddot(grad,grad)*relax;

           T = T2<Quantum>(t,false);
           compress(T,mpsxx::Left,0);
           compress(T,mpsxx::Right,0);

           eTA = exp(T,hf,2,D);
           normalize(eTA);

           E = inprod(mpsxx::Left,eTA,qc,eTA);

           cout << iter << "\t" << convergence << "\t" << E << endl;

           ++iter;

        }


     }

  template void gradient<Quantum>(const DArray<4> &,const MPO<Quantum> &,double E,const MPS<Quantum> &tccd,const MPS<Quantum> &wccd,const T2_2_mpo &list,DArray<4> &grad,bool merged);
  template void solve<Quantum>(DArray<4> &t,const MPO<Quantum> &qc,const MPS<Quantum> &hf,const std::vector<double> &,int,double);

}
