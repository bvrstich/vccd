#ifndef INPUT_H
#define INPUT_H

#include <iostream>
#include <cstdlib>

using std::ostream;

void read_oei(const char *,DArray<2> &,const std::vector<int> &);
void random_oei(DArray<2> &);
void read_tei(const char *,DArray<4> &,const std::vector<int> &);
void random_tei(DArray<4> &);

#endif
