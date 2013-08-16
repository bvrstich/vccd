#ifndef SPINHAMILTONIAN_H
#define SPINHAMILTONIAN_H

#include <iostream>
#include <iomanip>

//some functions which initialize an MPO to a certian Hamiltonian

using namespace btas;
using namespace mps;

template<class Q>
void physical(Qshapes<Q> &qp);

template<class Q>
MPO<Q> creator(int L,int site,int spin);

template<class Q>
MPO<Q> annihilator(int L,int site,int spin);

template<class Q>
MPO<Q> n_loc(int L,int site);

template<class Q>
MPO<Q> N_tot(int L);

template<class Q>
MPO<Q> n_up_tot(int L);

template<class Q>
MPO<Q> n_down_tot(int L);

template<class Q>
MPO<Q> hubbard(int L,double U);

template<class Q>
MPO<Q> T1(const DArray<2> &);

#endif
