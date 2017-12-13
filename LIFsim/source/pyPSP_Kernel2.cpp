#include "PSP_Kernel2.h"

#include <boost/python.hpp>
using namespace boost::python;

void export_PSP_Kernel2(){
    class_<PSP_Kernel2>("PSP_Kernel",init<double,double>())
        .def_readwrite("taum",&PSP_Kernel2::tau)
        .def_readwrite("taus",&PSP_Kernel2::tau_s)
        .def_readwrite("T",&PSP_Kernel2::T)
        .def_readonly("V0",&PSP_Kernel2::Vo)
        .def_readonly("lntt",&PSP_Kernel2::lntt)
        .def_readonly("t_max_fac",&PSP_Kernel2::t_max_fac)
        .def_readonly("v_max_fac",&PSP_Kernel2::v_max_fac)
        .def("__call__",&PSP_Kernel2::operator ())   
    ;
}

