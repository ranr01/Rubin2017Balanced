#include "Spike2.h"
#include <iostream>

void printSpike(Spike2 * sp){
    std::cout<<sp->affarent<<" "<<sp->time<<" "<<sp->timeBlock<<"\n";           
}

#include "boost/python.hpp"
using namespace boost::python;


void export_Spike2()
{
     def("printSpike",&printSpike);
     class_<Spike2,Spike2 *>("Spike2")
        .def(init<int,double>())
        .def(self < self)
        .def_readwrite("time",&Spike2::time)
        .def_readwrite("aff",&Spike2::affarent)
        .def_readwrite("timeBlock",&Spike2::timeBlock)
        .def_readwrite("tau_m_exp",&Spike2::tau_exponent)
        .def_readwrite("tau_s_exp",&Spike2::tau_s_exponent)
    ;    
}
