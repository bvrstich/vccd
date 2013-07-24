#include <iostream>
#include <iomanip>

using std::cout;
using std::endl;

#include "include.h"

using namespace btas;

/**
 * empty constructor
 */
Tensor::Tensor() : QSDArray<3>(){ }

/**
 * constructor will construct a QSDArray<3> and fill it randomly
 * @param qshape The allowed quantumnumbers corresponding to the three different legs of the tensor
 * @param dshape The corresponding dimensions of the three different legs of the tensor
 */
Tensor::Tensor(const blitz::TinyVector<Qshapes,3> &qshape,const blitz::TinyVector<Dshapes,3> &dshape) : QSDArray<3> (Quantum::zero(),qshape,dshape,Tools::rgen) { }

/**
 * copy constructor
 * @param tensor_c input Tensor to be copied
 */
Tensor::Tensor(const Tensor &tensor_c) : QSDArray<3> (tensor_c) { }

/**
 * destructor
 */
Tensor::~Tensor() {}

/**
 * fill the tensor randomly
 */
void Tensor::fill_Random(){

   *this = Tools::rgen;

}

/**
 * Canonicalize the tensor to a right or left normalized tensor
 */
void Tensor::canonicalize(bool left){

}
