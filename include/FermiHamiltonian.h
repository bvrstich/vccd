#ifndef FERMIHAMILTONIAN_H
#define FERMIHAMILTONIAN_H

#include <iostream>
#include <iomanip>

class Ostate;

//some functions which initialize an MPO to a certian Hamiltonian

using namespace btas;
using namespace mpsxx;

template<class Q>
void physical(Qshapes<Q> &qp);

template<class Q>
MPO<Q> crea_up(int,int);

template<class Q>
MPO<Q> crea_down(int,int);

template<class Q>
MPO<Q> anni_up(int,int);

template<class Q>
MPO<Q> anni_down(int,int);

template<class Q>
MPO<Q> crea_up_crea_up(int,int,int,double);

template<class Q>
MPO<Q> crea_up_crea_down(int,int,int,double);

template<class Q>
MPO<Q> crea_down_crea_up(int,int,int,double);

template<class Q>
MPO<Q> crea_down_crea_down(int,int,int,double);

template<class Q>
MPO<Q> anni_up_anni_up(int,int,int,double);

template<class Q>
MPO<Q> anni_down_anni_up(int,int,int,double);

template<class Q>
MPO<Q> anni_up_anni_down(int,int,int,double);

template<class Q>
MPO<Q> anni_down_anni_down(int,int,int,double);

template<class Q>
MPO<Q> E(int,int,int,double);

template<class Q>
MPO<Q> E(int,int,int,int,int,double);

template<class Q>
MPO<Q> tpint(int,int,int,int,int,double);

template<class Q>
MPO<Q> T1(const DArray<2> &,bool);

template<class Q>
MPO<Q> T2(const DArray<4> &,bool);

template<class Q>
MPO<Q> qcham(const DArray<2> &,const DArray<4> &,bool);

template<class Q>
MPO<Q> one_body(const DArray<2> &,bool);

template<class Q>
void get_merged_index(const Qshapes<Q> &,Qshapes<Q> &,std::vector< std::vector<int> > &,std::vector< std::vector<int> > &);

template<class Q>
int is_in(const Q &qn,const Qshapes<Q> &qlist);

void insert_id(QSDArray<4> &,int row,int column,double value = 1.0);
void insert_zero(QSDArray<4> &,int row,int column);
void insert_sign(QSDArray<4> &,int row,int column,double value = 1.0);
void insert_crea_up(QSDArray<4> &,int row,int column,double);
void insert_crea_up_s(QSDArray<4> &,int row,int column,double);
void insert_crea_down(QSDArray<4> &,int row,int column,double);
void insert_crea_down_s(QSDArray<4> &,int row,int column,double);
void insert_anni_up(QSDArray<4> &,int row,int column,double);
void insert_anni_up_s(QSDArray<4> &,int row,int column,double);
void insert_anni_down(QSDArray<4> &,int row,int column,double);
void insert_anni_down_s(QSDArray<4> &,int row,int column,double);
void insert_crea_up_crea_down(QSDArray<4> &,int row,int column,double);
void insert_crea_up_anni_up(QSDArray<4> &,int row,int column,double);
void insert_crea_up_anni_down(QSDArray<4> &,int row,int column,double);
void insert_crea_down_anni_up(QSDArray<4> &,int row,int column,double);
void insert_crea_down_anni_down(QSDArray<4> &,int row,int column,double);
void insert_anni_up_anni_down(QSDArray<4> &,int row,int column,double);
void insert_triple_crea_up_first(QSDArray<4> &,int row,int column,double,double);
void insert_triple_crea_up_last(QSDArray<4> &,int row,int column,double);
void insert_triple_crea_down_first(QSDArray<4> &,int row,int column,double,double);
void insert_triple_crea_down_last(QSDArray<4> &,int row,int column,double);
void insert_triple_anni_down_first(QSDArray<4> &,int row,int column,double,double);
void insert_triple_anni_down_last(QSDArray<4> &,int row,int column,double);
void insert_triple_anni_up_first(QSDArray<4> &,int row,int column,double,double);
void insert_triple_anni_up_last(QSDArray<4> &,int row,int column,double);
void insert_local(QSDArray<4> &,int row,int column,double,double);
void insert_pair(QSDArray<4> &,int row,int column,const std::vector<double> &);
void insert_pair_s(QSDArray<4> &,int row,int column,const std::vector<double> &);
void insert_local_ob(QSDArray<4> &,int row,int column,double,double);


void fill_mp2(DArray<4> &T,const DArray<4> &V,const std::vector<double> &ener);

#endif
