#ifndef SPIKE2_H
#define SPIKE2_H

#include <ostream>

class Spike2
{
 public:
  double time;
  long double tau_exponent, tau_s_exponent;
  int affarent, spike_index,timeBlock;
  Spike2() {}
  
  Spike2(int aff, double t, long double t_ex, long double ts_ex,int tb=0)
    : time(t), affarent(aff), tau_exponent(t_ex) , tau_s_exponent(ts_ex), timeBlock(tb) {}
  
  Spike2(int aff, double t)
        :time(t), affarent(aff),timeBlock(0){}
  
  friend std::ostream& operator<<(std::ostream& o , Spike2 const & sp);
  friend int operator<(Spike2 const & s1 , Spike2 const & s2);  
  friend int operator==(Spike2 const & s1 , Spike2 const & s2);
};





  
#endif
