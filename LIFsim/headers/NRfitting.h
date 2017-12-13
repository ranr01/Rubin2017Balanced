#ifndef NRFITTING_H
#define NRFITTING_H

#include <vector>

namespace NR{
  class model_function
    {
    public:
      std::vector<double> a;
      int ma;
      virtual double operator()(double x)=0;
      virtual std::vector<double> grad(double x)=0;
   
    };
  
  double chi_sqr(std::vector<double>& x,std::vector<double>& y,std::vector<double>& grad,std::vector< std::vector<double> >& hess,model_function & f);
  void mrqmin2D(std::vector<double>& x,std::vector<double>& y, model_function& f);

}
#endif
