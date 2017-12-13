#include <deque>
//#include <iostream>
#include <numeric>
#include "SpikingTempotron.h"
#include "Spike2.h"
#include "Tempotron2.h"
#include "NR.h"
//#include "InputLayer2.h"

using namespace std;

std::deque<MaxPoint> & SpikingTempotron::getVoltageTrace(){
  static long double lnD, lnD_s, t_max, V_cur, V_next;
    static std::deque< Spike2 >::const_iterator spikesEnd;
    static int dtb;
    
    KminusTh f(this);
    crossings_.clear();
    

    D_td_ = deque<long double>(spikes_->number_of_tds);
    D_s_td_ = deque<long double>(spikes_->number_of_tds);
    V_td_ = deque<long double>(spikes_->number_of_tds);

    V_t_.clear();
    
#ifdef DEBUG
    int stage;
#endif
    spikesEnd = spikes_->end();

    currentSpike_ = spikes_->begin();
    nextSpike_ = currentSpike_;
    ++nextSpike_;

     timeBlock_ = 0;
    expTimeBlock_taum_ = std::exp((long double)(-spikes_->timeBlockSize/K_->tau));
    expTimeBlock_taus_ = std::exp((long double)(-spikes_->timeBlockSize/K_->tau_s));
    TB_BeginTime_ = 0.0;
    reset_D_ = 0.0;
    D=0.0;
    D_s=0.0;

    try {
        while (nextSpike_ != spikesEnd) {
             if (currentSpike_->timeBlock > timeBlock_){
                    D *= std::pow(expTimeBlock_taum_,currentSpike_->timeBlock-timeBlock_);
                    reset_D_ *= std::pow(expTimeBlock_taum_,currentSpike_->timeBlock-timeBlock_);
                    D_s *= std::pow(expTimeBlock_taus_,currentSpike_->timeBlock-timeBlock_);
                    timeBlock_ = currentSpike_->timeBlock;
                    TB_BeginTime_ = timeBlock_ * spikes_->timeBlockSize;
                }
             
            D += w_[ currentSpike_->affarent ] * currentSpike_->tau_exponent;
            D_s += w_[ currentSpike_->affarent ] * currentSpike_->tau_s_exponent;
            V_t_.push_back(MaxPoint(TB_BeginTime_+currentSpike_->time,K_->Vo * 
                                                       (D / currentSpike_->tau_exponent - D_s / currentSpike_->tau_s_exponent)));
           
            if ((D > 0.0) && (D_s > 0.0)) {
                lnD = std::log(D);
                lnD_s = std::log(D_s);
                t_max = K_->t_max_fac * (lnD_s - lnD + K_->lntt);
                if (t_max > currentSpike_->time) {
                    dtb = nextSpike_->timeBlock - timeBlock_;
                    if (t_max <= nextSpike_->time+ dtb*spikes_->timeBlockSize) {
                        V_max_ = std::exp(K_->v_max_fac * (K_->tau * lnD - K_->tau_s * lnD_s));
                        if (V_max_ > threshold_) {
#ifdef DEBUG
                            stage = 1;
#endif
                            crossings_.push_back(TB_BeginTime_+NR::Math::find_x(currentSpike_->time, t_max, f));
                            V_t_.push_back(MaxPoint(crossings_.back(), threshold_));
                            reset_V_trace(crossings_.back());
                        }else
                                V_t_.push_back(MaxPoint(TB_BeginTime_+t_max, V_max_));
                    } else {
                        V_next = K_->Vo * (D / nextSpike_->tau_exponent * std::pow(expTimeBlock_taum_,dtb) -
                                    D_s / nextSpike_->tau_s_exponent * std::pow(expTimeBlock_taus_,dtb));
                        if (V_next > threshold_) {
#ifdef DEBUG
                            stage = 2;
#endif
                            crossings_.push_back(TB_BeginTime_+NR::Math::find_x(currentSpike_->time, nextSpike_->time+ dtb*spikes_->timeBlockSize, f));
                            V_t_.push_back(MaxPoint(crossings_.back(), threshold_));
                            reset_V_trace(crossings_.back());
                        }
                    }
                }
            }
            ++currentSpike_;
            ++nextSpike_;
        }

        if (currentSpike_->timeBlock > timeBlock_){
                    D *= std::pow(expTimeBlock_taum_,currentSpike_->timeBlock-timeBlock_);
                    reset_D_ *= std::pow(expTimeBlock_taum_,currentSpike_->timeBlock-timeBlock_);
                    D_s *= std::pow(expTimeBlock_taus_,currentSpike_->timeBlock-timeBlock_);
                    timeBlock_ = currentSpike_->timeBlock;
                    TB_BeginTime_ = timeBlock_ * spikes_->timeBlockSize;
                }
             
        D += w_[ currentSpike_->affarent ] * currentSpike_->tau_exponent;
        D_s += w_[ currentSpike_->affarent ] * currentSpike_->tau_s_exponent;
        V_t_.push_back(MaxPoint(TB_BeginTime_+currentSpike_->time, K_->Vo * (D / currentSpike_->tau_exponent - D_s / currentSpike_->tau_s_exponent)));
        
        if ((D > 0.0) && (D_s > 0.0)) {
            lnD = std::log(D);
            lnD_s = std::log(D_s);
            t_max = K_->t_max_fac * (lnD_s - lnD + K_->lntt);
            if (t_max > currentSpike_->time) {
                V_max_ = std::exp(K_->v_max_fac * (K_->tau * lnD - K_->tau_s * lnD_s));
                if (V_max_ > threshold_) {
#ifdef DEBUG
                    stage = 3;
#endif
                    crossings_.push_back(TB_BeginTime_+NR::Math::find_x(currentSpike_->time, t_max, f));
                    V_t_.push_back(MaxPoint(crossings_.back(), threshold_));
                    
                    SpikeTrain tmp;
                    tmp.push_back(Spike2(this->nSynapses_,1e6,1.,1.));
                    nextSpike_=tmp.begin();
                    reset_V_trace(crossings_.back());
                }
                V_t_.push_back(MaxPoint(TB_BeginTime_+t_max, V_max_));

            }
        }
    } catch (std::runtime_error & e) {
        cout << " Error in activate()\n";
#ifdef DEBUG
        cout << "stage=" << stage;
#endif
        cout << " S" << currentSpike_ - spikes_->begin() << " " << D << " " << D_s << "\n"
                << e.what() << "\n"
                << " t_max=" << t_max << " V_next=" << V_next << " V_max_=" << V_max_ << " thershold_=" << threshold_ << "\n";
        cout<< "curruntSpike: "<<*currentSpike_<<'\n'
             << "nextSpike: "<<*nextSpike_<<"\n"
             << "End=" << (nextSpike_==spikesEnd) << endl;
        throw;
    }

    ++currentSpike_;
    
    return V_t_;
    
}
void SpikingTempotron::output_trace_noTeacher(std::ostream & o) {
    static long double lnD, lnD_s, t_max, V_cur, V_next;
    static std::deque< Spike2 >::const_iterator spikesEnd;
    static int dtb;
    
    KminusTh f(this);
    crossings_.clear();

    D_td_ = deque<long double>(spikes_->number_of_tds);
    D_s_td_ = deque<long double>(spikes_->number_of_tds);
    V_td_ = deque<long double>(spikes_->number_of_tds);

#ifdef DEBUG
    int stage;
#endif
    spikesEnd = spikes_->end();

    currentSpike_ = spikes_->begin();
    nextSpike_ = currentSpike_;
    ++nextSpike_;

    reset_D_ = 0.0;
    D=0.0;
    D_s=0.0;

    try {
        while (nextSpike_ != spikesEnd) {
            D += w_[ currentSpike_->affarent ] * currentSpike_->tau_exponent;
            D_s += w_[ currentSpike_->affarent ] * currentSpike_->tau_s_exponent;
            if (currentSpike_->affarent==this->nSynapses_)
                o << MaxPoint(currentSpike_->time,K_->Vo * (D / currentSpike_->tau_exponent - D_s / currentSpike_->tau_s_exponent)) << " 3\n";
            else
                o << MaxPoint(currentSpike_->time, K_->Vo * (D / currentSpike_->tau_exponent - D_s / currentSpike_->tau_s_exponent)) << " 0\n";

            if ((D > 0.0) && (D_s > 0.0)) {
                lnD = std::log(D);
                lnD_s = std::log(D_s);
                t_max = K_->t_max_fac * (lnD_s - lnD + K_->lntt);
                if (t_max > currentSpike_->time) {
                    if (t_max <= nextSpike_->time) {
                        V_max_ = std::exp(K_->v_max_fac * (K_->tau * lnD - K_->tau_s * lnD_s));
                        if (V_max_ > threshold_) {
#ifdef DEBUG
                            stage = 1;
#endif
                            crossings_.push_back(NR::Math::find_x(currentSpike_->time, t_max, f));
                            o << MaxPoint(crossings_.back(), threshold_) << " 1\n";
                            reset_V_trace(crossings_.back(), o);
                        }
                        o << MaxPoint(t_max, V_max_) << " 2\n";
                    } else {
                        V_next = K_->Vo * (D / nextSpike_->tau_exponent - D_s / nextSpike_->tau_s_exponent);
                        if (V_next > threshold_) {
#ifdef DEBUG
                            stage = 2;
#endif
                            crossings_.push_back(NR::Math::find_x(currentSpike_->time, nextSpike_->time, f));
                            o << MaxPoint(crossings_.back(), threshold_) << " 1\n";
                            reset_V_trace(crossings_.back(), o);
                        }
                    }
                }
            }
            ++currentSpike_;
            ++nextSpike_;
        }

        D += w_[ currentSpike_->affarent ] * currentSpike_->tau_exponent;
        D_s += w_[ currentSpike_->affarent ] * currentSpike_->tau_s_exponent;
         if (currentSpike_->affarent==this->nSynapses_)
                o << MaxPoint(currentSpike_->time,K_->Vo * (D / currentSpike_->tau_exponent - D_s / currentSpike_->tau_s_exponent)) << " 3\n";
            else
                o << MaxPoint(currentSpike_->time, K_->Vo * (D / currentSpike_->tau_exponent - D_s / currentSpike_->tau_s_exponent)) << " 0\n";
        
        if ((D > 0.0) && (D_s > 0.0)) {
            lnD = std::log(D);
            lnD_s = std::log(D_s);
            t_max = K_->t_max_fac * (lnD_s - lnD + K_->lntt);
            if (t_max > currentSpike_->time) {
                V_max_ = std::exp(K_->v_max_fac * (K_->tau * lnD - K_->tau_s * lnD_s));
                if (V_max_ > threshold_) {
#ifdef DEBUG
                    stage = 3;
#endif
                    crossings_.push_back(NR::Math::find_x(currentSpike_->time, t_max, f));
                    o << MaxPoint(crossings_.back(), threshold_) << " 1\n";
                   reset_V_trace(crossings_.back(), o);
                }
                o << MaxPoint(t_max, V_max_) << " 2\n";

            }
        }
    } catch (std::runtime_error & e) {
        cout << " Error in activate()\n";
#ifdef DEBUG
        cout << "stage=" << stage;
#endif
        cout << " S" << currentSpike_ - spikes_->begin() << " " << D << " " << D_s << "\n"
                << e.what() << "\n"
                << " t_max=" << t_max << " V_next=" << V_next << " V_max_=" << V_max_ << " thershold_=" << threshold_ << "\n";
        cout<< "curruntSpike: "<<*currentSpike_<<'\n'
             << "nextSpike: "<<*nextSpike_<<"\n"
             << "End=" << (nextSpike_==spikesEnd) << endl;
        throw;
    }

    ++currentSpike_;
  
}
void SpikingTempotron::outputVoltageTrace(std::ostream & o) {
    static long double lnD, lnD_s, t_max, V_cur, V_next;
    static std::deque< Spike2 >::const_iterator spikesEnd;
    KminusTh f(this);
    crossings_.clear();
    D_td_ = deque<long double>(spikes_->number_of_tds);
    D_s_td_ = deque<long double>(spikes_->number_of_tds);
    V_td_ = deque<long double>(spikes_->number_of_tds);
    spikesEnd = spikes_->end();

    currentSpike_ = spikes_->begin();
    nextSpike_ = currentSpike_;
    ++nextSpike_;

    reset_D_ = 0.0;

    //cout<<"\ntd:"<<*tdSpike_<<"\ntd-eps:"<<*tdminusepsSpike_<<endl;
    for (td_count_ = 0; td_count_ < spikes_->number_of_tds; ++td_count_) {
        tdSpike_ = spikes_->tds[td_count_];
        tdminusepsSpike_ = spikes_->tds_minus_eps[td_count_];

        while (currentSpike_ != tdminusepsSpike_) {
            D += w_[ currentSpike_->affarent ] * currentSpike_->tau_exponent;
            D_s += w_[ currentSpike_->affarent ] * currentSpike_->tau_s_exponent;
            //cout<<"S"<<currentSpike_ - spikes_->begin() - 2<<" "<<D<<" "<<D_s<<" "<<flush;
            o << MaxPoint(currentSpike_->time, K_->Vo * (D / currentSpike_->tau_exponent - D_s / currentSpike_->tau_s_exponent)) << " 0\n";

            if ((D > 0.0) && (D_s > 0.0)) {
                lnD = std::log(D);
                lnD_s = std::log(D_s);
                t_max = K_->t_max_fac * (lnD_s - lnD + K_->lntt);
                if (t_max > currentSpike_->time) {
                    if (t_max <= nextSpike_->time) {
                        V_max_ = std::exp(K_->v_max_fac * (K_->tau * lnD - K_->tau_s * lnD_s));
                        if (V_max_ > threshold_) {
                            S_ = 1.0;
                            //std::cout<<"*1* "<<flush;
                            crossings_.push_back(NR::Math::find_x(currentSpike_->time, t_max, f));
                            o << MaxPoint(crossings_.back(), threshold_) << " 1\n";
                            reset_V_trace(crossings_.back(), o);
                        } else
                            o << MaxPoint(t_max, V_max_) << " 2\n";
                    } else {
                        V_next = K_->Vo * (D / nextSpike_->tau_exponent - D_s / nextSpike_->tau_s_exponent);
                        if (V_next > threshold_) {
                            S_ = 1.0;
                            //std::cout<<"*2* "<<flush;
                            crossings_.push_back(NR::Math::find_x(currentSpike_->time, nextSpike_->time, f));
                            o << MaxPoint(crossings_.back(), threshold_) << " 1\n";
                            reset_V_trace(crossings_.back(), o);
                        }
                    }
                }
            }
            ++currentSpike_;
            ++nextSpike_;
        }

        V_next = K_->Vo * (D / currentSpike_->tau_exponent - D_s / currentSpike_->tau_s_exponent);

        while (currentSpike_ != tdSpike_) {
            D += w_[ currentSpike_->affarent ] * currentSpike_->tau_exponent;
            D_s += w_[ currentSpike_->affarent ] * currentSpike_->tau_s_exponent;

            V_cur = V_next;

            o << MaxPoint(currentSpike_->time, V_cur) << " 0\n";

            V_next = K_->Vo * (D / nextSpike_->tau_exponent - D_s / nextSpike_->tau_s_exponent);

            if ((V_cur < threshold_) && ((D > 0.0) && (D_s > 0.0))) {
                lnD = std::log(D);
                lnD_s = std::log(D_s);
                t_max = K_->t_max_fac * (lnD_s - lnD + K_->lntt);
                if (t_max > currentSpike_->time) {
                    if (t_max <= nextSpike_->time) {
                        V_max_ = std::exp(K_->v_max_fac * (K_->tau * lnD - K_->tau_s * lnD_s));
                        if (V_max_ > threshold_) {
                            S_ = 1.0;
                            //std::cout<<"*3* "<<flush;
                            crossings_.push_back(NR::Math::find_x(currentSpike_->time, t_max, f));
                            o << MaxPoint(crossings_.back(), threshold_) << " 1\n";
                            o << MaxPoint(t_max, V_max_) << " 2\n";

                        }
                    } else {
                        if (V_next > threshold_) {
                            S_ = 1.0;
                            //std::cout<<"*4* "<<flush;
                            crossings_.push_back(NR::Math::find_x(currentSpike_->time, nextSpike_->time, f));
                            o << MaxPoint(crossings_.back(), threshold_) << " 1\n";
                        }
                    }
                }
            }
            ++currentSpike_;
            ++nextSpike_;
        }

        o << MaxPoint(currentSpike_->time, V_next) << " 3\n";

        teacher_reset_V(K_->Vo * (D / tdSpike_->tau_exponent - D_s / tdSpike_->tau_s_exponent));
    }
    while (nextSpike_ != spikesEnd) {
        D += w_[ currentSpike_->affarent ] * currentSpike_->tau_exponent;
        D_s += w_[ currentSpike_->affarent ] * currentSpike_->tau_s_exponent;

        if (currentSpike_ == tdSpike_)
            o << MaxPoint(currentSpike_->time, K_->Vo * (D / currentSpike_->tau_exponent - D_s / currentSpike_->tau_s_exponent)) << " 3\n";
        else
            o << MaxPoint(currentSpike_->time, K_->Vo * (D / currentSpike_->tau_exponent - D_s / currentSpike_->tau_s_exponent)) << " 0\n";

        if ((D > 0.0) && (D_s > 0.0)) {
            lnD = std::log(D);
            lnD_s = std::log(D_s);
            t_max = K_->t_max_fac * (lnD_s - lnD + K_->lntt);
            if (t_max > currentSpike_->time) {
                if (t_max <= nextSpike_->time) {
                    V_max_ = std::exp(K_->v_max_fac * (K_->tau * lnD - K_->tau_s * lnD_s));
                    if (V_max_ > threshold_) {
                        S_ = 1.0;
                        //std::cout<<"*5* "<<flush;
                        crossings_.push_back(NR::Math::find_x(currentSpike_->time, t_max, f));
                        o << MaxPoint(crossings_.back(), threshold_) << " 1\n";
                        reset_V_trace(crossings_.back(), o);
                    } else
                        o << MaxPoint(t_max, V_max_) << " 2\n";
                } else {
                    V_next = K_->Vo * (D / nextSpike_->tau_exponent - D_s / nextSpike_->tau_s_exponent);
                    if (V_next > threshold_) {
                        S_ = 1.0;
                        //std::cout<<"*6* "<<flush;
                        crossings_.push_back(NR::Math::find_x(currentSpike_->time, nextSpike_->time, f));
                        o << MaxPoint(crossings_.back(), threshold_) << " 1\n";
                        reset_V_trace(crossings_.back(), o);
                    }
                }
            }
        }
        ++currentSpike_;
        ++nextSpike_;
    }

    D += w_[ currentSpike_->affarent ] * currentSpike_->tau_exponent;
    D_s += w_[ currentSpike_->affarent ] * currentSpike_->tau_s_exponent;

    o << MaxPoint(currentSpike_->time, K_->Vo * (D / currentSpike_->tau_exponent - D_s / currentSpike_->tau_s_exponent)) << " 0\n";

    if ((D > 0.0) && (D_s > 0.0)) {
        lnD = std::log(D);
        lnD_s = std::log(D_s);
        t_max = K_->t_max_fac * (lnD_s - lnD + K_->lntt);
        if (t_max > currentSpike_->time) {
            V_max_ = std::exp(K_->v_max_fac * (K_->tau * lnD - K_->tau_s * lnD_s));
            if (V_max_ > threshold_) {
                S_ = 1.0;
                //std::cout<<"*7* "<<flush;
                crossings_.push_back(NR::Math::find_x(currentSpike_->time, t_max, f));
                o << MaxPoint(crossings_.back(), threshold_) << " 1\n";
            }
            o << MaxPoint(t_max, V_max_) << " 2\n";
        }
    }
    ++currentSpike_;
}

void SpikingTempotron::reset_V_trace(double t, ostream & o) {
    static long double lnD, lnD_s, t_max, V_next;
    static int dtb;
    
    KminusTh f(this);

    D -= (threshold_ / K_->Vo) * std::exp((long double) (t / K_->tau));

    o << MaxPoint(t, K_->Vo * (D * std::exp(-(long double) (t / K_->tau)) - D_s * std::exp(-(long double) (t / K_->tau_s)))) << " 1\n";

    if ((D > 0.0) && (D_s > 0.0)) {
        lnD = std::log(D);
        lnD_s = std::log(D_s);
        t_max = K_->t_max_fac * (lnD_s - lnD + K_->lntt);
        if (t_max > t) {
            dtb = nextSpike_->timeBlock - timeBlock_;
            if (t_max <= nextSpike_->time+ dtb*spikes_->timeBlockSize) {
                V_max_ = std::exp(K_->v_max_fac * (K_->tau * lnD - K_->tau_s * lnD_s));
                if (V_max_ > threshold_) {
                    t_max_ = t_max;
                    S_ = 1.0;
                    //std::cout<<"*1a* "<<flush;
                    crossings_.push_back(NR::Math::find_x(t, t_max, f));
                    o << MaxPoint(crossings_.back(), threshold_) << " 1\n";
                    reset_V_trace(crossings_.back(), o);
                } else
                    o << MaxPoint(t_max, V_max_) << " 2\n";
            } else {
                V_next = K_->Vo * (D / nextSpike_->tau_exponent * std::pow(expTimeBlock_taum_,dtb) -
                                    D_s / nextSpike_->tau_s_exponent * std::pow(expTimeBlock_taus_,dtb));
                if (V_next > threshold_) {
                    t_max_ = t_max;
                    S_ = 1.0;
                    //std::cout<<"*2a* "<<flush;
                    crossings_.push_back(NR::Math::find_x(t, nextSpike_->time+ dtb*spikes_->timeBlockSize, f));
                    o << MaxPoint(crossings_.back(), threshold_) << " 1\n";
                    reset_V_trace(crossings_.back(), o);
                }
            }
        }
    }

}
void SpikingTempotron::reset_V_trace(double t) {
    static long double lnD, lnD_s, t_max, V_next;
    static int dtb;

    KminusTh f(this);
    
    t -= TB_BeginTime_;
    
    D -= (fraction_reset*threshold_ / K_->Vo) * std::exp((long double) (t / K_->tau));

    V_t_.push_back(MaxPoint(TB_BeginTime_+t, K_->Vo * (D * std::exp((long double) (-t / K_->tau)) - D_s * std::exp((long double) (-t / K_->tau_s)))));

    if ((D > 0.0) && (D_s > 0.0)) {
        lnD = std::log(D);
        lnD_s = std::log(D_s);
        t_max = K_->t_max_fac * (lnD_s - lnD + K_->lntt);
        if (t_max > t) {
            dtb = nextSpike_->timeBlock - timeBlock_;
            if (t_max <= nextSpike_->time+ dtb*spikes_->timeBlockSize) {
                V_max_ = std::exp(K_->v_max_fac * (K_->tau * lnD - K_->tau_s * lnD_s));
                if (V_max_ > threshold_) {
                    t_max_ = t_max;
                    S_ = 1.0;
                    //std::cout<<"*1a* "<<flush;
                    crossings_.push_back(TB_BeginTime_+NR::Math::find_x(t, t_max, f));
                    V_t_.push_back(MaxPoint(crossings_.back(), threshold_));
                    reset_V_trace(crossings_.back());
                } else
                    V_t_.push_back(MaxPoint(TB_BeginTime_+t_max, V_max_));
            } else {
                V_next = K_->Vo * (D / nextSpike_->tau_exponent * std::pow(expTimeBlock_taum_,dtb) -
                                    D_s / nextSpike_->tau_s_exponent * std::pow(expTimeBlock_taus_,dtb));
                if (V_next > threshold_) {
                    t_max_ = t_max;
                    S_ = 1.0;
                    //std::cout<<"*2a* "<<flush;
                    crossings_.push_back(TB_BeginTime_+NR::Math::find_x(t, nextSpike_->time+ dtb*spikes_->timeBlockSize, f));
                    V_t_.push_back(MaxPoint(crossings_.back(), threshold_));
                    reset_V_trace(crossings_.back());
                }
            }
        }
    }

}
