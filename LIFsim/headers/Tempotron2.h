#ifndef TEMPOTRON2_H
#define TEMPOTRON2_H

#include <cstdlib>
#include <vector>
#include <deque>
#include <iterator>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <string>
#include "Spike2.h"
#include "PSP_Kernel2.h"
#include "SpikeTrain.h"
#include "MaxPoint.h"
#ifdef USE_PYTHON
        #include <boost/python.hpp>
#endif


class Tempotron2 
{
public:
  Tempotron2(int synapses, double lambda = 0.0001, double mu = 0.99, double threshold = 1.0);

  Tempotron2(Tempotron2 const & t);

  ~Tempotron2() {}

  virtual void activate();
    	   
  void learn_ltd();
  void learn_ltp();

   void restart();
   
   virtual void resetMomentum() { delta_w_gag_psik = std::vector<double>(nSynapses_+1,0.0); }
      
   void setSpikes(SpikeTrain  * s)
     { spikes_ = s; }
   void setPattern(SpikeTrain  * s)
     { spikes_ = s; }

   void setLmbda(double val) { lambda_ = val; }
   void setMu(double val) { mu_ = val; }
   void setKernel(PSP_Kernel2 * new_K) { K_ = new_K;
   #ifdef MAX_DT_EQ_T
    spike_at_T=Spike2(0, K_->T,std::exp(K_->T/K_->tau) , std::exp(K_->T/K_->tau_s) );
   #endif
   }
   PSP_Kernel2 * getKernel() {return K_;}
   void setThreshold(double val) { threshold_=val;}

   double getLmbda() { return lambda_ ; }
   double getMu() { return mu_ ; }  
   double getThreshold() { return threshold_; }
   double const & getVmax() const { return V_max_;}

   void rand_w(double std);       // Generates random w's normaly distributed ~N(0.0,std)
   std::vector<double>& w() {return w_;}
   double S() { return S_; }
   int getnSynapses() { return nSynapses_; }
   void loadVector(std::string vectorFileName, std::ostream& out = std::cout);
   void saveVector(std::string vectorFN , std::ostream& out = std::cout);
   void serializeVectors(std::ostream& o);
   void loadVectors(std::istream& i);

   double decisionTime() { return decision_time_; }
   double maxDTime() { return maxDTime_; }
   void maxDTime(double t) { maxDTime_ = t; };

//DEBUG: 
   void traceVoltage(double dt);

   std::deque< MaxPoint >& getVoltageTrace() { return voltageTrace_; } 

  //std::vector< double> & dw() { return dw_; }
   
   void set_w(std::vector<double> const & new_w);

#ifdef USE_PYTHON
   void pySet_w(boost::python::object const & o);
   std::vector<double>::const_iterator pyw_begin,pyw_end;
   PSP_Kernel2 pygetKernel(){return *K_;}
#endif
   void get_IR(std::deque< double > & spike_times_vec,std::deque< double > & spike_aff_vec,
               std::deque< long double > & D_vec,
               std::deque< long double > & Ds_vec,
               std::deque< double > & crossing_times_vec,
               double & t_max, double& V_max, bool is_debug=false);

   double get_V_max( bool is_debug=false);

   bool get_Crossing_time_plus(std::deque< double > & crossing_times_vec);
   void weightDecay(double lambda);
   
   /*void printState( std::ostream& o ) { o<<S_<<' '; }
   double tMax() { return t_max_; }
   double VMax() { return V_max_; }
   */
//END DEBUG 
   
   class KminusTh
   {
   protected:
       Tempotron2 * tmp;
   public:
       KminusTh(Tempotron2 * t) : tmp(t) {}
       double operator()(long double t) const;
   };

   bool restrict_weights_sign, allow_initial_inhibatory_weights;
  
   std::vector< double> w_, w_sign_;
   int max_aff_learn; 
   void find_max();
   
   

 protected:
   double lambda_,t_max_, V_max_ , mu_, S_ , threshold_, decision_time_ , maxDTime_;
   long double D, D_s ;
   int nSynapses_;
   SpikeTrain  * spikes_; // a vector to store the spike time of the input neurons
   std::deque< Spike2 >::const_iterator currentSpike_ , nextSpike_;
   std::deque< MaxPoint > maxPoints_ ;
   
   std::vector< double> delta_w_gag_psik; // Momentum vector
   PSP_Kernel2 * K_; //a Pointer to a PSP_Kernel
   #ifdef MAX_DT_EQ_T
   Spike2 spike_at_T;
   #endif

    void enable_weight_sign_restriction(){
        for (int i = 0; i < nSynapses_  ; i++) {
            w_sign_[i] = (w_[i]>0) ? 1.0 : -1.0;
        }


  }
 
  void learn( double sign) ;
  


//DEGUB:
   std::deque< MaxPoint > voltageTrace_ ;
//   std::vector< double > dw_;
//END DEBUG
#ifdef USE_PYTHON
   //friend void export_Tempotron();
#endif
};

#endif
