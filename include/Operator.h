#ifndef OPERATOR_H
#define OPERATOR_H

#include <iostream>
#include <cstdlib>
#include <vector>

using std::ostream;

/**
 * @author Brecht Verstichel
 * @date 24-08-2013
 * This class contains the different operators needed in the construction of the renormalized operators for the gradient
 */
class Operator {

   public:

      //constructor
      Operator();

      //copy constructor
      Operator(const Operator &);

      //destructor
      virtual ~Operator();

      static void init(const DArray<2> &,const DArray<4> &);

      static void clear();

      static void construct_id(int site,int row,int col,double val);

      static void construct_s(int site,int row,int col,double val);

      static void construct_cu(int site,int row,int col,double val);

      static void construct_cus(int site,int row,int col,double val);

      static void construct_cd(int site,int row,int col,double val);

      static void construct_cds(int site,int row,int col,double val);

      static void construct_au(int site,int row,int col,double val);

      static void construct_aus(int site,int row,int col,double val);

      static void construct_ad(int site,int row,int col,double val);

      static void construct_ads(int site,int row,int col,double val);

      static void construct_cucd(int site,int row,int col,double val);

      static void construct_cuau(int site,int row,int col,double val);

      static void construct_cuad(int site,int row,int col,double val);

      static void construct_cdau(int site,int row,int col,double val);

      static void construct_cdad(int site,int row,int col,double val);

      static void construct_adau(int site,int row,int col,double val);

      static void construct_tcuf(int site,int row,int col,double tval,double Vval);

      static void construct_tcul(int site,int row,int col,double Vval);

      static void construct_tcdf(int site,int row,int col,double tval,double Vval);

      static void construct_tcdl(int site,int row,int col,double Vval);

      static void construct_tauf(int site,int row,int col,double tval,double Vval);

      static void construct_taul(int site,int row,int col,double Vval);

      static void construct_tadf(int site,int row,int col,double tval,double Vval);

      static void construct_tadl(int site,int row,int col,double Vval);

      static void construct_local(int site,int row,int col,double tval,double Vval);

      static void construct_pair_s(int site,int row,int column,const std::vector<double> &val);

      static void construct_pair(int site,int row,int column,const std::vector<double> &val);

      static void print();

      static const QSDArray<2> &gop(int site,int row,int col);

      static bool gsparse(int site,int row,int col);

   private:

      //!physical dimension of the operators:
      static Qshapes<Quantum> qp;

      //!vector containing all the 4x4 operators on every site for incoming and outgoing.
      static std::vector< QSDArray<2> ** > op;

      //!vector containing the sparsity information
      static std::vector< bool* > sparse;

      //!number of terms on row and columns of MPO
      static std::vector< std::vector<int> > dim;

      //!number of sites
      static int L;

};

#endif
