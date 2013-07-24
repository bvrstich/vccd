#include <iostream>
#include <fstream>
#include <cmath>

#include "include.h"

using std::cout;
using std::endl;

/**
 * simple random number generator
 */
double Tools::rgen() { 
   
   return 2.0*(static_cast<double>(rand())/RAND_MAX) - 1.0; 
   
}
