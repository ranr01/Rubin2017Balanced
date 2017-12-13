#include <cstdlib>
#include <iostream>
#include "Spike2.h"

std::ostream& operator<<(std::ostream& o , Spike2 const & sp)
{ 
  o<<"Spike: t("<<sp.affarent<<")="<<sp.time<<" timeBlock="<<sp.timeBlock<<" "<<sp.tau_exponent<<" "<<sp.tau_s_exponent; 
  return o; 
}

int operator<(Spike2 const & s1 , Spike2 const & s2) { 
    return 
        (s1.timeBlock==s2.timeBlock) ? (s1.time < s2.time) : s1.timeBlock < s2.timeBlock; 
}
int operator==(Spike2 const & s1 , Spike2 const & s2) { return (s1.timeBlock==s2.timeBlock)&&(s1.time == s2.time); }
