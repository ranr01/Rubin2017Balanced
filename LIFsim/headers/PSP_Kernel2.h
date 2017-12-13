#ifndef PSP_KERNEL2_H
#define PSP_KERNEL2_H

#include <cstdlib>
#include <cmath>

class PSP_Kernel2
{
 public:
  long double tau, tau_s , T, v_max_fac;
  long double Vo, lntt, t_max_fac,  t_extr, Vo_inv;
  PSP_Kernel2(double t, double ts,double TT=1.0) : tau(t), tau_s(ts), T(TT)
  {
    t_max_fac = ( tau*tau_s )/( tau-tau_s );
    v_max_fac = 1.0 / ( tau-tau_s );
    Vo = v_max_fac * std::pow( tau , tau*v_max_fac ) * std::pow( tau_s , -tau_s*v_max_fac ) ;
    Vo_inv = 1.0/Vo;
    lntt = std::log( tau / tau_s );
    t_extr = t_max_fac * lntt;
  }

  double operator()(double t) { return (t>0) ? Vo*(std::exp(-t/tau) - std::exp(-t/tau_s)) : 0; }

};

#endif 
