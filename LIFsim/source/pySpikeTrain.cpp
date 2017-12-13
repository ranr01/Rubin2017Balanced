#include "SpikeTrain.h"
#include "RanUtilsRandom.h"
#include <stdexcept>
#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>


using namespace boost::python;

std::vector<double> * SpikeTrain::getTds(){
    std::vector<double> * a=new std::vector<double>(number_of_tds);
    for (int i=0 ; i<number_of_tds; ++i){
        (*a)[i]=tds[i]->timeBlock*timeBlockSize + tds[i]->time;
    }
    return a;
}
std::vector<double> * SpikeTrain::getTdsMinusEps(){
        std::vector<double> * a=new std::vector<double>(number_of_tds);
        for (int i=0 ; i<number_of_tds; ++i){
                (*a)[i]=tds_minus_eps[i]->timeBlock*timeBlockSize + tds_minus_eps[i]->time;
    }
    return a;
}

SpikeTrain * generatePoissonSpikeTrian(boost::python::object const & mean_rates,
                                       double T,int N,
                                       int N_noise){
    /* 
     * Generates a spike train with N Poisson spike times with mean rate mean_rates[i]
     * For Duration T. In addition adds a single spike at a random time from each 
     * noise afferent with N_noise noise afferents 
     */
    
    static RU::PoissonProccess p;
    static RU::Random rnd;
    static double t;
    SpikeTrain * st_ptr = new SpikeTrain();
    boost::python::stl_input_iterator<double> x_i(mean_rates),x_end;
    int i = 0;
    while (x_i!=x_end){
        // generate poisson rates with mean rate x_i
        p.setTau(1./(*x_i));
        t=0.0;
        while (t<T){
            t += p.NextEvent();
            st_ptr->add_spike(i,t);
        }
        st_ptr->pop_back();
        ++i;
        ++x_i;
    }
    if (i!=N)
        throw std::runtime_error("Number of elements in iterable object is not N");//assert i==N
    //generate noise spikes (1 spike per afferent.)
    for(i=N; i<N+N_noise; ++i){
        t = rnd.nextFloat(T);
        st_ptr->add_spike(i,t);
    }
    //sorting according to spike time   
    st_ptr->sort();
    
    return st_ptr;
}

SpikeTrain * generateNoisySpikeTrian(SpikeTrain * st, double sigma,
                                       double T,int N, int N_noise){
    static RU::Random rnd;
    NR::GausianDist g(0.0,sigma);
    SpikeTrain * st_ptr = new SpikeTrain();
    
    for (auto spike=st->begin(); spike!=st->end(); ++spike){
        st_ptr->add_spike(spike->affarent,spike->time+g());
    } 
    
    //generate noise spikes (1 spike per afferent.)
    for(int i=N; i<N+N_noise; ++i){
        st_ptr->add_spike(i,rnd.nextFloat(T));
    }
    //sorting according to spike time   
    st_ptr->sort();
    
    return st_ptr;
}


BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(setTdMarkers_OL,setTdMarkers,1,2)

void export_SpikeTrain(){
    class_<SpikeTrain>("SpikeTrain")
        .def(vector_indexing_suite<SpikeTrain>())
        .def("getTds",&SpikeTrain::getTds,return_value_policy<manage_new_object>())
        .def("getTdsMinusEps",&SpikeTrain::getTdsMinusEps,return_value_policy<manage_new_object>())
        .def("sort",&SpikeTrain::sort)
        .def("setTdMarkers",&SpikeTrain::setTdMarkers,setTdMarkers_OL())
        .def_readwrite("timeBlockSize",&SpikeTrain::timeBlockSize)
        .def("applyTimeBlock",&SpikeTrain::applyTimeBlock)
        .def("calcExponents",&SpikeTrain::calcExponents)
        .def("_cpp_add_spike",&SpikeTrain::add_spike)
    ;
    def("CPPgeneratePoissonSpikeTrian",generatePoissonSpikeTrian,return_value_policy<manage_new_object>());
    def("CPPgenerateNoisySpikeTrian",generateNoisySpikeTrian,return_value_policy<manage_new_object>());
}