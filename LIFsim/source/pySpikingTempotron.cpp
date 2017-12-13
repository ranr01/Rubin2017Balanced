#include "SpikingTempotron.h"
#include <stdexcept>
#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "vector_indexing_suite_no_contains.h"

using namespace boost::python;


void Tempotron2::pySet_w(boost::python::object const & o){
    boost::python::stl_input_iterator<double> begin(o),end;
    w_=std::vector<double>(begin,end);
    if (w_.size()!=nSynapses_)
           throw std::runtime_error("Number of element in iteratable object is not N");
    w_.push_back(0.0);
    w_.push_back(-0.01);
    pyw_begin=w_.begin();
    pyw_end=w_.end()-2;
    if (restrict_weights_sign)
        enable_weight_sign_restriction();
}



void export_Tempotron(){
    class_<MaxPoint>("__MaxPoint",init<double,double>())
        .def_readwrite("time",&MaxPoint::time)
        .def_readwrite("V",&MaxPoint::V)
    ;
    
    class_< std::deque<MaxPoint> >("__dequeMaxPoint")
        .def(vector_indexing_suite_no_contains< std::deque<MaxPoint>  >())
    ;
    
    class_<Tempotron2>("Tempotron",init<int, optional<double,double,double> >())
            .add_property("llambda",&Tempotron2::getLmbda,&Tempotron2::setLmbda)
            .add_property("mu",&Tempotron2::getMu,&Tempotron2::setMu)
            .add_property("threshold",&Tempotron2::getThreshold,&Tempotron2::setThreshold)
            .add_property("nSynapses",&Tempotron2::getnSynapses)
            .def("activate",&Tempotron2::activate)
            .def("rand_w",&Tempotron2::rand_w)
            .def("setPattern",&Tempotron2::setPattern)
            .add_property("K",&Tempotron2::pygetKernel,&Tempotron2::setKernel)
            .add_property("_w",range(&Tempotron2::pyw_begin,&Tempotron2::pyw_end),&Tempotron2::pySet_w)
            .def("restart",&Tempotron2::restart)
            .def("resetMomentum",&Tempotron2::resetMomentum)
            .def("S",&Tempotron2::S)
            .def_readwrite("_cpp_w",&Tempotron2::w_)
            .def_readwrite( "restrict_weights_sign",&Tempotron2::restrict_weights_sign)
            .def_readwrite("allow_initial_inhibatory_weights",&Tempotron2::allow_initial_inhibatory_weights)
            .def_readwrite("max_aff_learn",&Tempotron2::max_aff_learn)
            .def("_traceVoltage",&Tempotron2::traceVoltage)
            .def("_getVoltageTrace",&Tempotron2::getVoltageTrace,return_value_policy<reference_existing_object>())
            .def("__getw",&Tempotron2::w,return_value_policy<reference_existing_object>())
            .def("weightDecay",&Tempotron2::weightDecay)
            .def("get_V_max",&Tempotron2::get_V_max)
    ;
}

void export_SpikingTempotron(){

     class_<SpikingTempotron,bases<Tempotron2> >("SpikingTempotron",init<int, optional<double,double,double> >())
        .add_property("epsilon",&SpikingTempotron::getEpsilon,&SpikingTempotron::setEpsilon)
        .def("activate_no_teacher",&SpikingTempotron::activate_noTeacher)
        .def("crossings",&SpikingTempotron::crossings,return_value_policy<reference_existing_object>())
        .def("_getVoltageTrace",&SpikingTempotron::getVoltageTrace,return_value_policy<reference_existing_object>())
        .def_readwrite("FractionReset",&SpikingTempotron::fraction_reset)
#ifdef DEBUG
      .def_readonly("_stage",&SpikingTempotron::stage)
#endif
        ;

 
} 
