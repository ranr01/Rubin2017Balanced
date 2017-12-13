/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <vector>
#include <deque>

#include "NR.h"

void export_Spike2();
void export_PSP_Kernel2();
void export_SpikeTrain();
void export_Tempotron();
void export_SpikingTempotron();

void setRNGenSeed(long seed) { NR::Rnd.setSeed(seed); }
long getRNGenSeed() {return NR::Rnd.getSeed();}

using namespace boost::python;

template< template<class,class> class Container, class T >
void export_container(const char *  name){
    class_<Container<T,std::allocator<T> > >(name)
     .def(vector_indexing_suite< Container<T,std::allocator<T> > >())
     .def(init<int,T>())
     .def(init<int>())
    ;
}


BOOST_PYTHON_MODULE(_ST) {
    export_container<std::vector,double>("__vecdub");
    export_container<std::deque,double>("__dequedub");
    export_container<std::deque,long double>("__dequelongdouble");
    export_container<std::deque,short>("__dequeshort");
    export_container<std::vector, std::vector<double> >("__vecvecdouble");
    export_container<std::deque, std::vector<double> >("__dequevecdouble");
    export_container<std::deque, std::deque<double> >("__dequedequedouble");
     
    def("setRNGenSeed",&setRNGenSeed);
    def("getRNGenSeed",&getRNGenSeed);
    
    export_Spike2();
    export_PSP_Kernel2();
    export_SpikeTrain();
    export_Tempotron();
    export_SpikingTempotron();
}

