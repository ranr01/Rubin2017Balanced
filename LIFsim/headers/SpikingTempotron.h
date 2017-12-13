/* 
 * File:   SpikingTempotron.h
 * Author: ranr
 *
 * Created on April 29, 2010, 10:10 AM
 */

#ifndef _SPIKINGTEMPOTRON_H
#define	_SPIKINGTEMPOTRON_H

#include "SpikeTrain.h"
#include "Tempotron2.h"
#include <deque>
#include <vector>
#include <iostream>
//#include <boost/numeric/ublas/symmetric.hpp>

class SpikingTempotron : public Tempotron2 {
public:

    SpikingTempotron(int synapses, double lambda = 0.0001, double mu = 0.99, double threshold = 1.0)
    : Tempotron2(synapses, lambda, mu, threshold), dw_(synapses + 1, 0.0), dtheta_gag_(0.0), 
            multiple_spikes(0), expTimeBlock_taum_(1.), expTimeBlock_taus_(1.), timeBlock_(0), TB_BeginTime_ (0.),
    fraction_reset(1.) { 
    }

    ~SpikingTempotron() {
    }

    virtual void activate();

    void activate_noTeacher();
    virtual std::deque<MaxPoint> & getVoltageTrace();
    void outputVoltageTrace(std::ostream & o);
    void output_trace_noTeacher(std::ostream & o);

    virtual int learn();

    virtual void resetMomentum() {
        delta_w_gag_psik = std::vector<double>(nSynapses_+1,0.0);
        dtheta_gag_ = 0.0;
    }

    void setEpsilon(double eps) {
        epsilon_ = eps;
    }

    double getEpsilon() {
        return epsilon_;
    }

    std::vector< double > const & dw() {
        return dw_;
    }

    std::deque< double > const & crossings(){
        return crossings_;
    }

     int multiple_spikes;
     
     long double fraction_reset;
#ifdef DEBUG
     int stage;
#endif
protected:
    double epsilon_, dtheta_gag_;
    std::deque<double> crossings_;
    long double reset_D_;
    std::deque< long double > D_td_, D_s_td_, V_td_;
    std::deque< Spike2 >::const_iterator tdSpike_, tdminusepsSpike_;
    std::vector< double > dw_;
    int td_count_;
    std::deque<MaxPoint> V_t_;
    long double  expTimeBlock_taum_,expTimeBlock_taus_;
    int timeBlock_;
    double TB_BeginTime_;
    
    inline void changeTimeBlock(){
     
                D *= std::pow(expTimeBlock_taum_, currentSpike_->timeBlock - timeBlock_);
                reset_D_ *= std::pow(expTimeBlock_taum_, currentSpike_->timeBlock - timeBlock_);
                D_s *= std::pow(expTimeBlock_taus_, currentSpike_->timeBlock - timeBlock_);
                timeBlock_ = currentSpike_->timeBlock;
                TB_BeginTime_ = timeBlock_ * spikes_->timeBlockSize;
            
    }
    
    void reset_V(double t);
    void reset_V_last_spike(double t);
    void reset_V_trace(double t, std::ostream & o);
    virtual void reset_V_trace(double t);
    virtual void teacher_reset_V(long double Vtd);
    virtual bool td_err(int count=0);

};

#endif	/* _SPIKINGTEMPOTRON_H */

