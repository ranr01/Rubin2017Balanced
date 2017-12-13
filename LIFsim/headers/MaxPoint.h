/* 
 * File:   MaxPoint.h
 * Author: ranr
 *
 * Created on April 19, 2012, 11:05 PM
 */

#ifndef MAXPOINT_H
#define	MAXPOINT_H
#include <iostream>

class MaxPoint {
 public:
  double time, V;
  MaxPoint( double t, double v) : time(t), V(v) {}
  MaxPoint(){}
  friend std::ostream & operator<<(std::ostream & o , MaxPoint const & p);
  
#ifdef USE_PYTHON
 // bool operator==(MaxPoint const & other){
 //     return ( (time==other.time)&&(V==other.V) );
 // }
#endif
};


#endif	/* MAXPOINT_H */

