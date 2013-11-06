#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

using std::endl;
using std::ostream;
using std::ofstream;
using std::ifstream;

#include "include.h"

QSDArray<2> Operator::id;
QSDArray<2> Operator::s;
QSDArray<2> Operator::cu;
QSDArray<2> Operator::cd;
QSDArray<2> Operator::cus;
QSDArray<2> Operator::cds;
QSDArray<2> Operator::au;
QSDArray<2> Operator::ad;
QSDArray<2> Operator::aus;
QSDArray<2> Operator::ads;
QSDArray<2> Operator::cucd;
QSDArray<2> Operator::cuau;
QSDArray<2> Operator::cuad;
QSDArray<2> Operator::cdau;
QSDArray<2> Operator::cdad;
QSDArray<2> Operator::auad;
QSDArray<2> Operator::tcu;
QSDArray<2> Operator::tcd;
QSDArray<2> Operator::tad;
QSDArray<2> Operator::tau;
QSDArray<2> Operator::lob;
QSDArray<2> Operator::ltb;

//standard constructor
Operator::Operator(){ }

Operator::Operator(const Operator &){ }

Operator::~Operator(){ } 

/**
 * initialize all the operators
 */
void Operator::init(){ 
  
   Qshapes<Quantum> qp;
   physical(qp);

   Quantum qt;
   DArray<2> Ip(1,1);
   DArray<2> Im(1,1);

   Ip = 1.0;
   Im = -1.0;

   //first id
   qt = Quantum::zero();

   id.resize(qt,make_array(qp,-qp));

   id.insert(shape(0,0),Ip);
   id.insert(shape(1,1),Ip);
   id.insert(shape(2,2),Ip);
   id.insert(shape(3,3),Ip);

   //sign
   s.resize(qt,make_array(qp,-qp));

   Ip = 1.0;

   s.insert(shape(0,0),Ip);
   s.insert(shape(1,1),Im);
   s.insert(shape(2,2),Im);
   s.insert(shape(3,3),Ip);

   //create up
   qt = Quantum(1,0);

   cu.resize(qt,make_array(qp,-qp));

   cu.insert(shape(1,0),Ip);
   cu.insert(shape(3,2),Ip);

   //create up sign
   cus.resize(qt,make_array(qp,-qp));

   cus.insert(shape(1,0),Ip);
   cus.insert(shape(3,2),Im);

   //create down
   qt = Quantum(0,1);

   cd.resize(qt,make_array(qp,-qp));

   cd.insert(shape(2,0),Ip);
   cd.insert(shape(3,1),Ip);

   //create down sign
   cds.resize(qt,make_array(qp,-qp));

   cds.insert(shape(2,0),Ip);
   cds.insert(shape(3,1),Im);

   //annihilator of up spin
   qt = Quantum(-1,0);

   au.resize(qt,make_array(qp,-qp));

   au.insert(shape(0,1),Ip);
   au.insert(shape(2,3),Ip);

   //annihilator of up spin sign
   aus.resize(qt,make_array(qp,-qp));

   aus.insert(shape(0,1),Ip);
   aus.insert(shape(2,3),Im);

   //annihilator of down spin
   qt = Quantum(0,-1);

   ad.resize(qt,make_array(qp,-qp));

   ad.insert(shape(0,2),Ip);
   ad.insert(shape(1,3),Ip);

   //annihilator of down spin with sign
   ads.resize(qt,make_array(qp,-qp));

   ads.insert(shape(0,2),Ip);
   ads.insert(shape(1,3),Im);

   //creator of an up down pair on site
   qt = Quantum(1,1);

   cucd.resize(qt,make_array(qp,-qp));

   cucd.insert(shape(3,0),Ip);

   //create up annihilate up
   qt = Quantum::zero();

   cuau.resize(qt,make_array(qp,-qp));

   cuau.insert(shape(1,1),Ip);
   cuau.insert(shape(3,3),Ip);

   //create up annihilate down
   qt = Quantum(1,-1);

   cuad.resize(qt,make_array(qp,-qp));

   cuad.insert(shape(1,2),Ip);

   //create down annihilate up
   qt = Quantum(-1,1);

   cdau.resize(qt,make_array(qp,-qp));

   cdau.insert(shape(2,1),Ip);

   //create down annihilate down --> n_down
   qt = Quantum::zero();

   cdad.resize(qt,make_array(qp,-qp));

   cdad.insert(shape(2,2),Ip);
   cdad.insert(shape(3,3),Ip);

   //insert annihilate up annihilate down: extra minus sign! because it should be annihilate down annihilate up
   qt = Quantum(-1,-1);

   auad.resize(qt,make_array(qp,-qp));

   auad.insert(shape(0,3),Im);

   //complementary operator for a_up: behaves as a creator of up
   qt = Quantum(1,0);

   tcu.resize(qt,make_array(qp,-qp));

   tcu.insert(shape(3,2),Ip);

   //complementary operator for a_down: behaves as a creator of down 
   qt = Quantum(0,1);

   tcd.resize(qt,make_array(qp,-qp));

   tcd.insert(shape(3,1),Ip);

   //complementary operator for a^+_down: behaves as an annihilator of down 
   qt = Quantum(0,-1);

   tad.resize(qt,make_array(qp,-qp));

   tad.insert(shape(1,3),Ip);

   //complementary operator for a^+_up: behaves as an annihilator of up
   qt = Quantum(-1,0);

   tau.resize(qt,make_array(qp,-qp));

   tau.insert(shape(2,3),Ip);

   //local two body term
   qt = Quantum::zero();

   ltb.resize(qt,make_array(qp,-qp));

   ltb.insert(shape(3,3),Ip);

   //local one body term
   lob.resize(qt,make_array(qp,-qp));

   lob.insert(shape(1,1),Ip);
   lob.insert(shape(2,2),Ip);

   Ip = 2.0;

   lob.insert(shape(3,3),Ip);

}

/**
 * print all the operators, for testing purposes
 */
void Operator::print(){

   cout << "Identity" << endl;
   cout << endl;
   cout << Operator::id << endl;
   cout << endl;

   cout << "Sign" << endl;
   cout << endl;
   cout << Operator::s << endl;
   cout << endl;

   cout << "create up" << endl;
   cout << endl;
   cout << Operator::cu << endl;
   cout << endl;

   cout << "create up sign" << endl;
   cout << endl;
   cout << Operator::cus << endl;
   cout << endl;

   cout << "create down" << endl;
   cout << endl;
   cout << Operator::cd << endl;
   cout << endl;

   cout << "create down sign" << endl;
   cout << endl;
   cout << Operator::cds << endl;
   cout << endl;

   cout << "annihilate up" << endl;
   cout << endl;
   cout << Operator::au << endl;
   cout << endl;

   cout << "annihilate up sign" << endl;
   cout << endl;
   cout << Operator::aus << endl;
   cout << endl;

   cout << "annihilate down" << endl;
   cout << endl;
   cout << Operator::ad << endl;
   cout << endl;

   cout << "annihilate down isgn" << endl;
   cout << endl;
   cout << Operator::ads << endl;
   cout << endl;

   cout << "cu cd" << endl;
   cout << endl;
   cout << Operator::cucd << endl;
   cout << endl;

   cout << "cu au" << endl;
   cout << endl;
   cout << Operator::cuau << endl;
   cout << endl;

   cout << "cu ad" << endl;
   cout << endl;
   cout << Operator::cuad << endl;
   cout << endl;

   cout << "cd au" << endl;
   cout << endl;
   cout << Operator::cdau << endl;
   cout << endl;

   cout << "cd ad" << endl;
   cout << endl;
   cout << Operator::cdad << endl;
   cout << endl;

   cout << "au ad" << endl;
   cout << endl;
   cout << Operator::auad << endl;
   cout << endl;

   cout << "tcu" << endl;
   cout << endl;
   cout << Operator::tcu << endl;
   cout << endl;

   cout << "tcd" << endl;
   cout << endl;
   cout << Operator::tcd << endl;
   cout << endl;

   cout << "tad" << endl;
   cout << endl;
   cout << Operator::tad << endl;
   cout << endl;

   cout << "tau" << endl;
   cout << endl;
   cout << Operator::tau << endl;
   cout << endl;

   cout << "local two body" << endl;
   cout << endl;
   cout << Operator::ltb << endl;
   cout << endl;

   cout << "local one body" << endl;
   cout << endl;
   cout << Operator::lob << endl;
   cout << endl;

}
