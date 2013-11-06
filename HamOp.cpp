#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

using std::endl;
using std::ostream;
using std::ofstream;
using std::ifstream;

#include "include.h"

std::vector< Qshapes<Quantum> > HamOp::qo;

std::vector< std::vector< Ostate > > HamOp::ostates;

int HamOp::L;

/**
 * standard constructor
 */
HamOp::HamOp(){ }

/**
 * copy constructor
 */
HamOp::HamOp(const HamOp &ho_c){ }

/**
 * destructor
 */
HamOp::~HamOp(){ }

/**
 * initialize the static lists containing the incoming and outgoing states and quantumnumbers
 */
void HamOp::init(int L_in){

   L = L_in;

   qo.resize(L - 1);
   ostates.resize(L - 1);

   Ostate state;

   //outgoing states identity
   state.push_id();
   ostates[0].push_back(state);
   state.clear();

   qo[0].push_back(Quantum::zero());//I

   //singles

   //a^+_up
   state.push_crea_up(0);
   ostates[0].push_back(state);
   state.clear();

   qo[0].push_back(Quantum(-1,0));

   //a^+_down 
   state.push_crea_down(0);
   ostates[0].push_back(state);
   state.clear();

   qo[0].push_back(Quantum(0,-1));

   //a_up 
   state.push_anni_up(0);
   ostates[0].push_back(state);
   state.clear();

   qo[0].push_back(Quantum(1,0));

   //a_down
   state.push_anni_down(0);
   ostates[0].push_back(state);
   state.clear();

   qo[0].push_back(Quantum(0,1));

   //doubles:

   //a^+_up a^+_down
   state.push_crea_down(0);
   state.push_crea_up(0);
   ostates[0].push_back(state);
   state.clear();

   qo[0].push_back(Quantum(-1,-1));

   //doubles: a^+_up a_up
   state.push_anni_up(0);
   state.push_crea_up(0);
   ostates[0].push_back(state);
   state.clear();

   qo[0].push_back(Quantum::zero());

   //doubles: a^+_up a_down
   state.push_anni_down(0);
   state.push_crea_up(0);
   ostates[0].push_back(state);
   state.clear();

   qo[0].push_back(Quantum(-1,1));

   //doubles: a^+_down a_up
   state.push_anni_up(0);
   state.push_crea_down(0);
   ostates[0].push_back(state);
   state.clear();

   qo[0].push_back(Quantum(1,-1));

   //doubles: a^+_down a_down
   state.push_anni_down(0);
   state.push_crea_down(0);
   ostates[0].push_back(state);
   state.clear();

   qo[0].push_back(Quantum::zero());

   //doubles: a_up a_down
   state.push_anni_up(0);
   state.push_anni_down(0);
   ostates[0].push_back(state);
   state.clear();

   qo[0].push_back(Quantum(1,1));

   //complementary operators: triples: they have the state signature of the operator they are going to, but the opposite quantumnumber
   for(int j = 1;j < L;++j){

      state.push_crea_up(j);
      ostates[0].push_back(state);
      state.clear();

      qo[0].push_back(Quantum(1,0));

      state.push_crea_down(j);
      ostates[0].push_back(state);
      state.clear();

      qo[0].push_back(Quantum(0,1));

      state.push_anni_up(j);
      ostates[0].push_back(state);
      state.clear();

      qo[0].push_back(Quantum(-1,0));

      state.push_anni_down(j);
      ostates[0].push_back(state);
      state.clear();

      qo[0].push_back(Quantum(0,-1));

   }

   //finally the local term:
   state.push_id();
   ostates[0].push_back(state);
   state.clear();

   qo[0].push_back(Quantum::zero());

   //middle tensors
   for(int i = 1;i < L/2;++i){

      //identity
      state.push_id();
      ostates[i].push_back(state);
      state.clear();

      qo[i].push_back(Quantum::zero());//I

      //singles

      //a^+_up
      state.push_crea_up(i);
      ostates[i].push_back(state);
      state.clear();

      qo[i].push_back(Quantum(-1,0));

      //a^+_down 
      state.push_crea_down(i);
      ostates[i].push_back(state);
      state.clear();

      qo[i].push_back(Quantum(0,-1));

      //a_up 
      state.push_anni_up(i);
      ostates[i].push_back(state);
      state.clear();

      qo[i].push_back(Quantum(1,0));

      //a_down
      state.push_anni_down(i);
      ostates[i].push_back(state);
      state.clear();

      qo[i].push_back(Quantum(0,1));

      //copy the singles from previous tensor
      int row = 1;

      while(ostates[i - 1][row].size() == 1){

         ostates[i].push_back(ostates[i - 1][row]);
         qo[i].push_back(qo[i - 1][row]);
         ++row;

      }

      //doubles:

      //a^+_up a^+_down
      state.push_crea_down(i);
      state.push_crea_up(i);
      ostates[i].push_back(state);
      state.clear();

      qo[i].push_back(Quantum(-1,-1));

      //doubles: a^+_up a_up
      state.push_anni_up(i);
      state.push_crea_up(i);
      ostates[i].push_back(state);
      state.clear();

      qo[i].push_back(Quantum::zero());

      //doubles: a^+_up a_down
      state.push_anni_down(i);
      state.push_crea_up(i);
      ostates[i].push_back(state);
      state.clear();

      qo[i].push_back(Quantum(-1,1));

      //doubles: a^+_down a_up
      state.push_anni_up(i);
      state.push_crea_down(i);
      ostates[i].push_back(state);
      state.clear();

      qo[i].push_back(Quantum(1,-1));

      //doubles: a^+_down a_down
      state.push_anni_down(i);
      state.push_crea_down(i);
      ostates[i].push_back(state);
      state.clear();

      qo[i].push_back(Quantum::zero());

      //doubles: a_down a_up 
      state.push_anni_up(i);
      state.push_anni_down(i);
      ostates[i].push_back(state);
      state.clear();

      qo[i].push_back(Quantum(1,1));

      //make new doubles by adding on single f^'s and f's
      row = 1;

      while(ostates[i - 1][row].size() == 1){//create up

         state.push_crea_up(i);
         state.insert(state.end(),ostates[i - 1][row].begin(),ostates[i - 1][row].end());
         ostates[i].push_back(state);
         state.clear();

         Quantum tmp = qo[i - 1][row];
         tmp.crea_up();
         qo[i].push_back(tmp);

         ++row;

      }

      row = 1;

      while(ostates[i - 1][row].size() == 1){//create down

         state.push_crea_down(i);
         state.insert(state.end(),ostates[i - 1][row].begin(),ostates[i - 1][row].end());
         ostates[i].push_back(state);
         state.clear();

         Quantum tmp = qo[i - 1][row];
         tmp.crea_down();
         qo[i].push_back(tmp);

         ++row;

      }

      row = 1;

      while(ostates[i - 1][row].size() == 1){//anni up

         state.push_anni_up(i);
         state.insert(state.end(),ostates[i - 1][row].begin(),ostates[i - 1][row].end());
         ostates[i].push_back(state);
         state.clear();

         Quantum tmp = qo[i - 1][row];
         tmp.anni_up();
         qo[i].push_back(tmp);

         ++row;

      }

      row = 1;

      while(ostates[i - 1][row].size() == 1){//anni down 

         state.push_anni_down(i);
         state.insert(state.end(),ostates[i - 1][row].begin(),ostates[i - 1][row].end());
         ostates[i].push_back(state);
         state.clear();

         Quantum tmp = qo[i - 1][row];
         tmp.anni_down();
         qo[i].push_back(tmp);

         ++row;

      }

      //copy the doubles from previous tensor
      while(ostates[i - 1][row].size() == 2){

         ostates[i].push_back(ostates[i - 1][row]);
         qo[i].push_back(qo[i - 1][row]);
         ++row;

      }

      //complementary operators: triples
      for(int j = i + 1;j < L;++j){

         state.push_crea_up(j);
         ostates[i].push_back(state);
         state.clear();

         qo[i].push_back(Quantum(1,0));

         state.push_crea_down(j);
         ostates[i].push_back(state);
         state.clear();

         qo[i].push_back(Quantum(0,1));

         state.push_anni_up(j);
         ostates[i].push_back(state);
         state.clear();

         qo[i].push_back(Quantum(-1,0));

         state.push_anni_down(j);
         ostates[i].push_back(state);
         state.clear();

         qo[i].push_back(Quantum(0,-1));

      }

      //finally the local term:
      state.push_id();
      ostates[i].push_back(state);
      state.clear();

      qo[i].push_back(Quantum::zero());

   }

   //now switch around incoming and outgoing pairs to reduce the dimension of the MPO

   //first row is closed
   state.push_id();
   ostates[L/2].push_back(state);
   state.clear();

   qo[L/2].push_back(Quantum::zero());//I

   //outgoing singles
   for(int j = L/2 + 1;j < L;++j){

      state.push_crea_up(j);
      ostates[L/2].push_back(state);
      state.clear();

      qo[L/2].push_back(Quantum(1,0));

      state.push_crea_down(j);
      ostates[L/2].push_back(state);
      state.clear();

      qo[L/2].push_back(Quantum(0,1));

      state.push_anni_up(j);
      ostates[L/2].push_back(state);
      state.clear();

      qo[L/2].push_back(Quantum(-1,0));

      state.push_anni_down(j);
      ostates[L/2].push_back(state);
      state.clear();

      qo[L/2].push_back(Quantum(0,-1));

   }

   //outgoing pairs: inverse quantum numbers!
   Ostate istate;

   for(int i = L/2 + 1;i < L;++i){

      //first crea up
      istate.push_crea_up(i);

      //first on-site terms

      //a^+_up a^+_down
      state = istate;
      state.push_crea_down(i);
      ostates[L/2].push_back(state);
      state.clear();

      qo[L/2].push_back(Quantum(1,1));

      //a^+_up a_up
      state = istate;
      state.push_anni_up(i);
      ostates[L/2].push_back(state);
      state.clear();

      qo[L/2].push_back(Quantum::zero());

      //a^+_up a_down
      state = istate;
      state.push_anni_down(i);
      ostates[L/2].push_back(state);
      state.clear();

      qo[L/2].push_back(Quantum(1,-1));

      //two-site terms
      for(int j = i + 1;j < L;++j){

         //a^+_up a^+_up
         state = istate;
         state.push_crea_up(j);
         ostates[L/2].push_back(state);
         state.clear();

         qo[L/2].push_back(Quantum(2,0));

         //a^+_up a^+_down
         state = istate;
         state.push_crea_down(j);
         ostates[L/2].push_back(state);
         state.clear();

         qo[L/2].push_back(Quantum(1,1));

         //a^+_up a_down
         state = istate;
         state.push_anni_down(j);
         ostates[L/2].push_back(state);
         state.clear();

         qo[L/2].push_back(Quantum(1,-1));

         //a^+_up a_up
         state = istate;
         state.push_anni_up(j);
         ostates[L/2].push_back(state);
         state.clear();

         qo[L/2].push_back(Quantum(0,0));

      }

      istate.clear();

      //first crea down
      istate.push_crea_down(i);

      //a^+_down a_up
      state = istate;
      state.push_anni_up(i);
      ostates[L/2].push_back(state);
      state.clear();

      qo[L/2].push_back(Quantum(-1,1));

      //a^+_down a_down
      state = istate;
      state.push_anni_down(i);
      ostates[L/2].push_back(state);
      state.clear();

      qo[L/2].push_back(Quantum::zero());

      for(int j = i + 1;j < L;++j){

         //a^+_down a^+_up
         state = istate;
         state.push_crea_up(j);
         ostates[L/2].push_back(state);
         state.clear();

         qo[L/2].push_back(Quantum(1,1));

         //a^+_down a^+_down
         state = istate;
         state.push_crea_down(j);
         ostates[L/2].push_back(state);
         state.clear();

         qo[L/2].push_back(Quantum(0,2));

         //a^+_down a_down
         state = istate;
         state.push_anni_down(j);
         ostates[L/2].push_back(state);
         state.clear();

         qo[L/2].push_back(Quantum(0,0));

         //a^+_down a_up
         state = istate;
         state.push_anni_up(j);
         ostates[L/2].push_back(state);
         state.clear();

         qo[L/2].push_back(Quantum(-1,1));

      }

      istate.clear();

      //first anni down
      istate.push_anni_down(i);

      //a_down a_up 
      state = istate;
      state.push_anni_up(i);
      ostates[L/2].push_back(state);
      state.clear();

      qo[L/2].push_back(Quantum(-1,-1));

      for(int j = i + 1;j < L;++j){

         //a_down a^+_up
         state = istate;
         state.push_crea_up(j);
         ostates[L/2].push_back(state);
         state.clear();

         qo[L/2].push_back(Quantum(1,-1));

         //a_down a^+_down
         state = istate;
         state.push_crea_down(j);
         ostates[L/2].push_back(state);
         state.clear();

         qo[L/2].push_back(Quantum(0,0));

         //a_down a_down
         state = istate;
         state.push_anni_down(j);
         ostates[L/2].push_back(state);
         state.clear();

         qo[L/2].push_back(Quantum(0,-2));

         //a_down a_up
         state = istate;
         state.push_anni_up(j);
         ostates[L/2].push_back(state);
         state.clear();

         qo[L/2].push_back(Quantum(-1,-1));

      }
      
      istate.clear();

      //first is anni up 
      istate.push_anni_up(i);

      for(int j = i + 1;j < L;++j){

         //a_up a^+_up
         state = istate;
         state.push_crea_up(j);
         ostates[L/2].push_back(state);
         state.clear();

         qo[L/2].push_back(Quantum(0,0));

         //a_up a^+_down
         state = istate;
         state.push_crea_down(j);
         ostates[L/2].push_back(state);
         state.clear();

         qo[L/2].push_back(Quantum(-1,1));

         //a_up a_down
         state = istate;
         state.push_anni_down(j);
         ostates[L/2].push_back(state);
         state.clear();

         qo[L/2].push_back(Quantum(-1,-1));

         //a_up a_up
         state = istate;
         state.push_anni_up(j);
         ostates[L/2].push_back(state);
         state.clear();

         qo[L/2].push_back(Quantum(-2,0));

      }

      istate.clear();

   }

   //incoming singles
   for(int i = 0;i < L/2 + 1;++i){

      //a^+_up
      state.push_crea_up(i);
      ostates[L/2].push_back(state);
      state.clear();

      qo[L/2].push_back(Quantum(-1,0));

      //a^+_down
      state.push_crea_down(i);
      ostates[L/2].push_back(state);
      state.clear();

      qo[L/2].push_back(Quantum(0,-1));

      //a_down
      state.push_anni_down(i);
      ostates[L/2].push_back(state);
      state.clear();

      qo[L/2].push_back(Quantum(0,1));

      //a_up
      state.push_anni_up(i);
      ostates[L/2].push_back(state);
      state.clear();

      qo[L/2].push_back(Quantum(1,0));

   }

   //finally last column the identity
   state.push_id();
   ostates[L/2].push_back(state);
   state.clear();

   qo[L/2].push_back(Quantum::zero());

   //the rest of the blocks: L/2 + 1 --> L - 1
   for(int i = L/2 + 1;i < L - 1;++i){

      //first row is closed
      state.push_id();
      ostates[i].push_back(state);
      state.clear();

      qo[i].push_back(Quantum::zero());//I

      //outgoing singles
      for(int j = i + 1;j < L;++j){

         state.push_crea_up(j);
         ostates[i].push_back(state);
         state.clear();

         qo[i].push_back(Quantum(1,0));

         state.push_crea_down(j);
         ostates[i].push_back(state);
         state.clear();

         qo[i].push_back(Quantum(0,1));

         state.push_anni_up(j);
         ostates[i].push_back(state);
         state.clear();

         qo[i].push_back(Quantum(-1,0));

         state.push_anni_down(j);
         ostates[i].push_back(state);
         state.clear();

         qo[i].push_back(Quantum(0,-1));

      }

      //outgoing pairs: inverse quantum numbers!
      Ostate istate;

      for(int j = i + 1;j < L;++j){

         //first crea up
         istate.push_crea_up(j);

         //first on-site terms

         //a^+_up a^+_down
         state = istate;
         state.push_crea_down(j);
         ostates[i].push_back(state);
         state.clear();

         qo[i].push_back(Quantum(1,1));

         //a^+_up a_up
         state = istate;
         state.push_anni_up(j);
         ostates[i].push_back(state);
         state.clear();

         qo[i].push_back(Quantum::zero());

         //a^+_up a_down
         state = istate;
         state.push_anni_down(j);
         ostates[i].push_back(state);
         state.clear();

         qo[i].push_back(Quantum(1,-1));

         //two-site terms
         for(int k = j + 1;k < L;++k){

            //a^+_up a^+_up
            state = istate;
            state.push_crea_up(k);
            ostates[i].push_back(state);
            state.clear();

            qo[i].push_back(Quantum(2,0));

            //a^+_up a^+_down
            state = istate;
            state.push_crea_down(k);
            ostates[i].push_back(state);
            state.clear();

            qo[i].push_back(Quantum(1,1));

            //a^+_up a_down
            state = istate;
            state.push_anni_down(k);
            ostates[i].push_back(state);
            state.clear();

            qo[i].push_back(Quantum(1,-1));

            //a^+_up a_up
            state = istate;
            state.push_anni_up(k);
            ostates[i].push_back(state);
            state.clear();

            qo[i].push_back(Quantum(0,0));

         }

         istate.clear();

         //first crea down
         istate.push_crea_down(j);

         //a^+_down a_up
         state = istate;
         state.push_anni_up(j);
         ostates[i].push_back(state);
         state.clear();

         qo[i].push_back(Quantum(-1,1));

         //a^+_down a_down
         state = istate;
         state.push_anni_down(j);
         ostates[i].push_back(state);
         state.clear();

         qo[i].push_back(Quantum::zero());

         for(int k = j + 1;k < L;++k){

            //a^+_down a^+_up
            state = istate;
            state.push_crea_up(k);
            ostates[i].push_back(state);
            state.clear();

            qo[i].push_back(Quantum(1,1));

            //a^+_down a^+_down
            state = istate;
            state.push_crea_down(k);
            ostates[i].push_back(state);
            state.clear();

            qo[i].push_back(Quantum(0,2));

            //a^+_down a_down
            state = istate;
            state.push_anni_down(k);
            ostates[i].push_back(state);
            state.clear();

            qo[i].push_back(Quantum(0,0));

            //a^+_down a_up
            state = istate;
            state.push_anni_up(k);
            ostates[i].push_back(state);
            state.clear();

            qo[i].push_back(Quantum(-1,1));

         }

         istate.clear();

         //first anni down
         istate.push_anni_down(j);

         //a_down a_up 
         state = istate;
         state.push_anni_up(j);
         ostates[i].push_back(state);
         state.clear();

         qo[i].push_back(Quantum(-1,-1));

         for(int k = j + 1;k < L;++k){

            //a_down a^+_up
            state = istate;
            state.push_crea_up(k);
            ostates[i].push_back(state);
            state.clear();

            qo[i].push_back(Quantum(1,-1));

            //a_down a^+_down
            state = istate;
            state.push_crea_down(k);
            ostates[i].push_back(state);
            state.clear();

            qo[i].push_back(Quantum(0,0));

            //a_down a_down
            state = istate;
            state.push_anni_down(k);
            ostates[i].push_back(state);
            state.clear();

            qo[i].push_back(Quantum(0,-2));

            //a_down a_up
            state = istate;
            state.push_anni_up(k);
            ostates[i].push_back(state);
            state.clear();

            qo[i].push_back(Quantum(-1,-1));

         }

         istate.clear();

         //first is anni up 
         istate.push_anni_up(j);

         for(int k = j + 1;k < L;++k){

            //a_up a^+_up
            state = istate;
            state.push_crea_up(k);
            ostates[i].push_back(state);
            state.clear();

            qo[i].push_back(Quantum(0,0));

            //a_up a^+_down
            state = istate;
            state.push_crea_down(k);
            ostates[i].push_back(state);
            state.clear();

            qo[i].push_back(Quantum(-1,1));

            //a_up a_down
            state = istate;
            state.push_anni_down(k);
            ostates[i].push_back(state);
            state.clear();

            qo[i].push_back(Quantum(-1,-1));

            //a_up a_up
            state = istate;
            state.push_anni_up(k);
            ostates[i].push_back(state);
            state.clear();

            qo[i].push_back(Quantum(-2,0));

         }

         istate.clear();

      }

      //incoming singles
      for(int j = 0;j <= i;++j){

         //a^+_up
         state.push_crea_up(j);
         ostates[i].push_back(state);
         state.clear();

         qo[i].push_back(Quantum(-1,0));

         //a^+_down
         state.push_crea_down(j);
         ostates[i].push_back(state);
         state.clear();

         qo[i].push_back(Quantum(0,-1));

         //a_down
         state.push_anni_down(j);
         ostates[i].push_back(state);
         state.clear();

         qo[i].push_back(Quantum(0,1));

         //a_up
         state.push_anni_up(j);
         ostates[i].push_back(state);
         state.clear();

         qo[i].push_back(Quantum(1,0));

      }

      //finally last column the identity
      state.push_id();
      ostates[i].push_back(state);
      state.clear();

      qo[i].push_back(Quantum::zero());

   }

}

/**
 * print the states at the different boundaries
 */
void HamOp::print_states(){

   for(int i = 0;i < L - 1;++i){

      cout << "right boundary on site " << i << endl;
      cout << endl;
      for(int j = 0;j < ostates[i].size();++j){

         cout << endl;
         cout << ostates[i][j] << qo[i][j] << endl;
         cout << endl;

      }

   }

}


