#include <iostream>
#include <cstdlib>
#include <vector>
#include <deque>
#include <cmath>
#include <fstream>
//#include <bits/stl_deque.h>
//#include "RanUtils.h"

#include "PSP_Kernel2.h"
#include "Tempotron2.h"
#include "NR.h"

using namespace std;
//using namespace RU;

ostream & operator<<(ostream & o, MaxPoint const & p) {
    o << p.time << ' ' << p.V;
    return o;
}

Tempotron2::Tempotron2(int synapses, double lambda, double mu, double threshold) :
w_(synapses + 2, 0.0),
w_sign_(synapses),
delta_w_gag_psik(synapses + 1, 0.0),
t_max_(0),
lambda_(lambda),
mu_(mu),
threshold_(threshold),
restrict_weights_sign(false),
allow_initial_inhibatory_weights(true),
max_aff_learn(synapses)
{
    S_ = 0.0;
    nSynapses_ = synapses;
    w_[synapses + 1] = -0.01;

}

Tempotron2::Tempotron2(Tempotron2 const & t) : t_max_(0) {
    nSynapses_ = t.nSynapses_;
    lambda_ = t.lambda_;
    mu_ = t.mu_;
    threshold_ = t.threshold_;
    w_ = t.w_;
    w_sign_ = t.w_sign_;
    delta_w_gag_psik = t.delta_w_gag_psik;
    S_ = 0.0;
    K_ = t.K_;
    restrict_weights_sign = t.restrict_weights_sign;
    allow_initial_inhibatory_weights = t.allow_initial_inhibatory_weights;
    restart();
}

bool Tempotron2::get_Crossing_time_plus(std::deque< double > & crossing_times_vec) {
    restart();

    static long double lnD, lnD_s, t_star, V_star, V, V_next;
    std::deque< Spike2 >::const_iterator spikesEnd = spikes_->end();
    static KminusTh f(this);
    currentSpike_ = spikes_->begin();
    nextSpike_ = spikes_->begin();
    ++nextSpike_;

    V_next = -0.01;

    while ((nextSpike_ != spikesEnd)) {
        D += w_[ currentSpike_->affarent ] * currentSpike_->tau_exponent;
        D_s += w_[ currentSpike_->affarent ] * currentSpike_->tau_s_exponent;

        V = V_next;
        V_next = K_->Vo * (D / nextSpike_->tau_exponent - D_s / nextSpike_->tau_s_exponent);

        if ((D > 0.0) && (D_s > 0.0)) {
            lnD = std::log(D);
            lnD_s = std::log(D_s);
            t_star = K_->t_max_fac * (lnD_s - lnD + K_->lntt);
            if ((t_star > currentSpike_->time) && (t_star <= nextSpike_->time)) {
                V_star = std::exp(K_->v_max_fac * (K_->tau * lnD - K_->tau_s * lnD_s));
                if (V > threshold_) {
                    if (V_next <= threshold_) {
                        crossing_times_vec.push_back(NR::Math::find_x(t_star, nextSpike_->time, f));
                    }
                } else {
                    if (V_star > threshold_) {
                        S_ = 1.0;
                        crossing_times_vec.push_back(NR::Math::find_x(currentSpike_->time, t_star, f));
                        if (V_next <= threshold_) {
                            crossing_times_vec.push_back(NR::Math::find_x(t_star, nextSpike_->time, f));
                        }
                    }
                }
            } else {
                if (((V > threshold_) && (V_next <= threshold_)) || ((V <= threshold_) && (V_next > threshold_))) {
                    S_ = 1.0;
                    crossing_times_vec.push_back(NR::Math::find_x(currentSpike_->time, nextSpike_->time, f));
                }
            }
        } else if ((V > threshold_) && (V_next <= threshold_)) {
            crossing_times_vec.push_back(NR::Math::find_x(currentSpike_->time, nextSpike_->time, f));
        }
        ++currentSpike_;
        ++nextSpike_;
    }

    D += w_[ currentSpike_->affarent ] * currentSpike_->tau_exponent;
    D_s += w_[ currentSpike_->affarent ] * currentSpike_->tau_s_exponent;

    V = V_next;

    if ((D > 0.0) && (D_s > 0.0)) {
        lnD = std::log(D);
        lnD_s = std::log(D_s);
        t_star = K_->t_max_fac * (lnD_s - lnD + K_->lntt);
        if (t_star > currentSpike_->time) {
            V_star = std::exp(K_->v_max_fac * (K_->tau * lnD - K_->tau_s * lnD_s));
            if (V > threshold_) {
                crossing_times_vec.push_back(NR::Math::find_x(t_star, t_star + 10 * K_->tau, f));
            } else {
                if (V_star > threshold_) {
                    S_ = 1.0;
                    crossing_times_vec.push_back(NR::Math::find_x(currentSpike_->time, t_star, f));
                    crossing_times_vec.push_back(NR::Math::find_x(t_star, t_star + 10 * K_->tau, f));
                }
            }
        } else if (V > threshold_) {
            crossing_times_vec.push_back(NR::Math::find_x(currentSpike_->time, currentSpike_->time + 10 * K_->tau, f));
        }

    } else if (V > threshold_) {
        crossing_times_vec.push_back(NR::Math::find_x(currentSpike_->time, currentSpike_->time + 10 * K_->tau, f));
    }

    return (S_ == 1.0);

}

void Tempotron2::get_IR(std::deque< double > & spike_times_vec, std::deque< double > & spike_aff_vec,
        std::deque< long double > & D_vec,
        std::deque< long double > & Ds_vec,
        std::deque< double > & crossing_times_vec, double & t_max, double& V_max, bool is_debug) {
    restart();

    static long double lnD, lnD_s, t_star, V_star, V, V_next;
    std::deque< Spike2 >::const_iterator spikesEnd = spikes_->end();
    currentSpike_ = spikes_->begin();
    nextSpike_ = spikes_->begin();
    ++nextSpike_;

    V_max = -0.01;
    t_max = 0.0;

    V_next = -0.01;

    while ((nextSpike_ != spikesEnd)) {
        D += w_[ currentSpike_->affarent ] * currentSpike_->tau_exponent;
        D_s += w_[ currentSpike_->affarent ] * currentSpike_->tau_s_exponent;

        V = V_next;
        //if (is_debug) std::cout<<"V("<<currentSpike_-spikes_->begin()<<")="<<V<<" ";
        V_next = K_->Vo * (D / nextSpike_->tau_exponent - D_s / nextSpike_->tau_s_exponent);

        spike_times_vec.push_back(currentSpike_->time);
        spike_aff_vec.push_back(w_[ currentSpike_->affarent ]);
        D_vec.push_back(D);
        Ds_vec.push_back(D_s);

        V_max = (V < V_max) ? V_max : V;
        t_max = (V < V_max) ? t_max : currentSpike_->time;

        if ((D > 0.0) && (D_s > 0.0)) {
            lnD = std::log(D);
            lnD_s = std::log(D_s);
            t_star = K_->t_max_fac * (lnD_s - lnD + K_->lntt);
            if ((t_star > currentSpike_->time) && (t_star <= nextSpike_->time)) {
                V_star = std::exp(K_->v_max_fac * (K_->tau * lnD - K_->tau_s * lnD_s));
                V_max = (V_star < V_max) ? V_max : V_star;
                t_max = (V_star < V_max) ? t_max : t_star;
                if (V > threshold_) {
                    if (V_next <= threshold_) {
                        //if (is_debug) std::cout<<"*1* ";
                        KminusTh f(this);
                        crossing_times_vec.push_back(NR::Math::find_x(t_star, nextSpike_->time, f));
                    }
                } else {
                    if (V_star > threshold_) {
                        //if (is_debug) std::cout<<"*2* "<<"V_star="<<V_star<<" V="<<V<<" V_next="<<V_next<<" ";
                        KminusTh f(this);
                        crossing_times_vec.push_back(NR::Math::find_x(currentSpike_->time, t_star, f));
                        if (V_next <= threshold_) {
                            //if (is_debug) std::cout<<"*3* ";
                            crossing_times_vec.push_back(NR::Math::find_x(t_star, nextSpike_->time, f));
                        }

                    }
                }
            } else {
                if (((V > threshold_) && (V_next <= threshold_)) || ((V <= threshold_) && (V_next > threshold_))) {
                    KminusTh f(this);
                    //if (is_debug) std::cout<<"*4* ";
                    crossing_times_vec.push_back(NR::Math::find_x(currentSpike_->time, nextSpike_->time, f));
                }
            }
        } else if ((V > threshold_) && (V_next <= threshold_)) {
            KminusTh f(this);
            //if (is_debug) std::cout<<"*5* ";
            crossing_times_vec.push_back(NR::Math::find_x(currentSpike_->time, nextSpike_->time, f));
        }
        ++currentSpike_;
        ++nextSpike_;
    }

    // last spike
    D += w_[ currentSpike_->affarent ] * currentSpike_->tau_exponent;
    D_s += w_[ currentSpike_->affarent ] * currentSpike_->tau_s_exponent;

    V = V_next;
    //if (is_debug) std::cout<<"V("<<currentSpike_-spikes_->begin()<<")="<<V<<" D="<<D<<" D_s="<<D_s<<" ";

    spike_times_vec.push_back(currentSpike_->time);
    D_vec.push_back(D);
    Ds_vec.push_back(D_s);

    V_max = (V < V_max) ? V_max : V;
    t_max = (V < V_max) ? t_max : currentSpike_->time;

    if ((D > 0.0) && (D_s > 0.0)) {
        lnD = std::log(D);
        lnD_s = std::log(D_s);
        t_star = K_->t_max_fac * (lnD_s - lnD + K_->lntt);
        //if (is_debug) std::cout<<"t="<<currentSpike_->time<<" t_star="<<t_star<<" ";
        if (t_star > currentSpike_->time) {
            V_star = std::exp(K_->v_max_fac * (K_->tau * lnD - K_->tau_s * lnD_s));
            V_max = (V_star < V_max) ? V_max : V_star;
            t_max = (V_star < V_max) ? t_max : t_star;
            if (V > threshold_) {
                KminusTh f(this);
                //if (is_debug) std::cout<<"*6* ";
                crossing_times_vec.push_back(NR::Math::find_x(t_star, t_star + 10 * K_->tau, f));
            } else {
                if (V_star > threshold_) {
                    KminusTh f(this);
                    //if (is_debug) std::cout<<"*7* ";
                    crossing_times_vec.push_back(NR::Math::find_x(currentSpike_->time, t_star, f));
                    //if (is_debug) std::cout<<"*8* ";
                    crossing_times_vec.push_back(NR::Math::find_x(t_star, t_star + 10 * K_->tau, f));
                }
            }
        } else if (V > threshold_) {
            KminusTh f(this);
            //if (is_debug) std::cout<<"*9* ";
            crossing_times_vec.push_back(NR::Math::find_x(currentSpike_->time, currentSpike_->time + 10 * K_->tau, f));
        }

    } else if (V > threshold_) {
        KminusTh f(this);
        //if (is_debug) std::cout<<"*10* ";
        crossing_times_vec.push_back(NR::Math::find_x(currentSpike_->time, currentSpike_->time + 10 * K_->tau, f));
    }

}

double Tempotron2::get_V_max( bool is_debug) {
    restart();

    static long double lnD, lnD_s, t_star, V_star, V, V_next;
    std::deque< Spike2 >::const_iterator spikesEnd = spikes_->end();
    currentSpike_ = spikes_->begin();
    nextSpike_ = spikes_->begin();
    ++nextSpike_;

    long double V_max = -0.01;
    double t_max = 0.0;

    V_next = -0.01;

    while ((nextSpike_ != spikesEnd)) {
        D += w_[ currentSpike_->affarent ] * currentSpike_->tau_exponent;
        D_s += w_[ currentSpike_->affarent ] * currentSpike_->tau_s_exponent;

        V = V_next;
        //if (is_debug) std::cout<<"V("<<currentSpike_-spikes_->begin()<<")="<<V<<" ";
        V_next = K_->Vo * (D / nextSpike_->tau_exponent - D_s / nextSpike_->tau_s_exponent);

        V_max = (V < V_max) ? V_max : V;
        t_max = (V < V_max) ? t_max : currentSpike_->time;

        if ((D > 0.0) && (D_s > 0.0)) {
            lnD = std::log(D);
            lnD_s = std::log(D_s);
            t_star = K_->t_max_fac * (lnD_s - lnD + K_->lntt);
            if ((t_star > currentSpike_->time) && (t_star <= nextSpike_->time)) {
                V_star = std::exp(K_->v_max_fac * (K_->tau * lnD - K_->tau_s * lnD_s));
                V_max = (V_star < V_max) ? V_max : V_star;
                t_max = (V_star < V_max) ? t_max : t_star;   
            } 
        }
        ++currentSpike_;
        ++nextSpike_;
    }

    // last spike
    D += w_[ currentSpike_->affarent ] * currentSpike_->tau_exponent;
    D_s += w_[ currentSpike_->affarent ] * currentSpike_->tau_s_exponent;

    V = V_next;
    //if (is_debug) std::cout<<"V("<<currentSpike_-spikes_->begin()<<")="<<V<<" D="<<D<<" D_s="<<D_s<<" ";

    V_max = (V < V_max) ? V_max : V;
    t_max = (V < V_max) ? t_max : currentSpike_->time;

    if ((D > 0.0) && (D_s > 0.0)) {
        lnD = std::log(D);
        lnD_s = std::log(D_s);
        t_star = K_->t_max_fac * (lnD_s - lnD + K_->lntt);
        //if (is_debug) std::cout<<"t="<<currentSpike_->time<<" t_star="<<t_star<<" ";
        if (t_star > currentSpike_->time) {
            V_star = std::exp(K_->v_max_fac * (K_->tau * lnD - K_->tau_s * lnD_s));
            V_max = (V_star < V_max) ? V_max : V_star;
            t_max = (V_star < V_max) ? t_max : t_star;
            
        }
    }
    
    return double(V_max);

}

void Tempotron2::activate() {
    static long double lnD, lnD_s, t_max, V_next;
    static std::deque< Spike2 >::const_iterator spikesEnd;

    spikesEnd = spikes_->end();
    currentSpike_ = spikes_->begin();
    nextSpike_ = spikes_->begin();
    ++nextSpike_;

    while ((nextSpike_ != spikesEnd) && (S_ != 1.0)) {
        //cout<<currentSpike_->affarent<<' '<<flush;
        D += w_[ currentSpike_->affarent ] * currentSpike_->tau_exponent;
        D_s += w_[ currentSpike_->affarent ] * currentSpike_->tau_s_exponent;

        if ((D > 0.0) && (D_s > 0.0)) {
            lnD = std::log(D);
            lnD_s = std::log(D_s);
            t_max = K_->t_max_fac * (lnD_s - lnD + K_->lntt);
            if (t_max > currentSpike_->time) {
                if (t_max <= nextSpike_->time) {
                    V_max_ = std::exp(K_->v_max_fac * (K_->tau * lnD - K_->tau_s * lnD_s));
                    if (V_max_ > threshold_) {
                        t_max_ = t_max;
                        S_ = 1.0;

#ifdef TEMPOTRON2_SPIKE_TIME
                        KminusTh f(this);
                        decision_time_ = NR::Math::find_x(currentSpike_->time, t_max, f);
#endif
                    } else maxPoints_.push_back(MaxPoint(t_max, V_max_));
                } else {
                    V_next = K_->Vo * (D / nextSpike_->tau_exponent - D_s / nextSpike_->tau_s_exponent);
                    if (V_next > threshold_) {
                        t_max_ = t_max;
                        S_ = 1.0;
#if defined (TEMPOTRON2_SPIKE_TIME) || defined (ENABLE_T_MAX)
                        KminusTh f(this);
                        decision_time_ = NR::Math::find_x(currentSpike_->time, nextSpike_->time, f);
#endif
                        V_max_ = std::exp( K_->v_max_fac * ( K_->tau*lnD - K_->tau_s*lnD_s ) );
                        
                    } else maxPoints_.push_back(MaxPoint(nextSpike_->time, V_next));
                }
            }
        }
        ++currentSpike_;
        ++nextSpike_;
    }

#ifdef MAX_DT_EQ_T
    if (S_ != 1.0) {
        D += w_[ currentSpike_->affarent ] * currentSpike_->tau_exponent;
        D_s += w_[ currentSpike_->affarent ] * currentSpike_->tau_s_exponent;

        if ((D > 0.0) && (D_s > 0.0)) {
            lnD = std::log(D);
            lnD_s = std::log(D_s);
            t_max = K_->t_max_fac * (lnD_s - lnD + K_->lntt);
            if (t_max > currentSpike_->time) {
                V_max_ = std::exp(K_->v_max_fac * (K_->tau * lnD - K_->tau_s * lnD_s));
                if (V_max_ > threshold_) {
                    if (t_max <= K_->T) {
                        t_max_ = t_max;
                        S_ = 1.0;
#if defined (TEMPOTRON2_SPIKE_TIME) || defined (ENABLE_T_MAX)
                        KminusTh f(this);
                        decision_time_ = NR::Math::find_x(currentSpike_->time, t_max, f);
#endif
                    } else {
                        V_next = K_->Vo * (D / spike_at_T.tau_exponent - D_s / spike_at_T.tau_s_exponent);
                        if (V_next > threshold_) {
                            t_max_ = t_max;
                            S_ = 1.0;
#if defined (TEMPOTRON2_SPIKE_TIME) || defined (ENABLE_T_MAX)
                            KminusTh f(this);
                            decision_time_ = NR::Math::find_x(currentSpike_->time, nextSpike_->time, f);
#endif
                            //DEBUG
                            //V_max_ = std::exp( K_->v_max_fac * ( K_->tau*lnD - K_->tau_s*lnD_s ) );
                            //END DEBUG
                        }

                    }
                } else maxPoints_.push_back(MaxPoint(t_max, V_max_));
            }
        }
        ++currentSpike_;
    }
#else
    if (S_ != 1.0) {
        D += w_[ currentSpike_->affarent ] * currentSpike_->tau_exponent;
        D_s += w_[ currentSpike_->affarent ] * currentSpike_->tau_s_exponent;

        if ((D > 0.0) && (D_s > 0.0)) {
            lnD = std::log(D);
            lnD_s = std::log(D_s);
            t_max = K_->t_max_fac * (lnD_s - lnD + K_->lntt);
            if (t_max > currentSpike_->time) {
                V_max_ = std::exp(K_->v_max_fac * (K_->tau * lnD - K_->tau_s * lnD_s));
                if (V_max_ > threshold_) {
                    t_max_ = t_max;
                    S_ = 1.0;
#if defined (TEMPOTRON2_SPIKE_TIME) || defined (ENABLE_T_MAX)
                    KminusTh f(this);
                    decision_time_ = NR::Math::find_x(currentSpike_->time, t_max, f);
#endif
                } else maxPoints_.push_back(MaxPoint(t_max, V_max_));
            }
        }
        ++currentSpike_;
    }
#endif
    if (S_ == 0.0)
        find_max();

}

void Tempotron2::find_max() {
    V_max_ = -0.01;
    t_max_ = 0.0;

#ifdef ENABLE_T_MAX
    for (int i = 0; (i < maxPoints_.size()) && (maxPoints_[i].time < maxDTime_); ++i)
#else
    for (int i = 0; i < maxPoints_.size(); ++i)
#endif
    {

        if (maxPoints_[i].V > V_max_) {
            V_max_ = maxPoints_[i].V;
            t_max_ = maxPoints_[i].time;
        }
    }

    if (V_max_ <= 0.0) {
        t_max_ = (spikes_->rbegin())->time + K_->t_extr;
#ifdef ENABLE_T_MAX
        if (t_max_ > maxDTime_)
            t_max_ = maxDTime_;
#endif
        currentSpike_ = spikes_->end();
    } else {
        currentSpike_ = std::lower_bound(spikes_->begin(), spikes_->end(), Spike2(0, t_max_, 0, 0));
        if (currentSpike_ == spikes_->begin()) {
            t_max_ = (spikes_->rbegin())->time + K_->t_extr;
            currentSpike_ = spikes_->end();
        }
    }


    //std::cout<<maxPoints_.size()<<' '<<"V("<<maxPoints_[0].time<<")="<<maxPoints_[0].V<<' '<<"V("<<maxPoints_[1].time<<")="<<maxPoints_[1].V<<' '<<t_max_<<' '<<*currentSpike_<<' '<< (currentSpike_==spikes_->begin()) <<' '<<(currentSpike_==spikes_->end())<<' '<<std::flush;
    //std::cout<<"2learn(1) "<<std::flush;

}

void Tempotron2::learn_ltd() {
    //std::cout<<"2learn(-1) "<<std::flush;
#ifdef ENABLE_T_MAX
    t_max_ = (t_max_ > maxDTime_) ? maxDTime_ : t_max_;
#endif

    learn(-1.0);
}

void Tempotron2::learn_ltp() {

    learn(1.0);
}

void Tempotron2::learn(double sign) {
    //removed from this code since no learning is needed
}

void Tempotron2::restart() {
    t_max_ = 0.0;
    S_ = 0.0;
    D = 0.0;
    D_s = 0.0;
    maxPoints_ = deque< MaxPoint > ();
    decision_time_ = -0.1;
}

void Tempotron2::set_w(std::vector<double> const & new_w) {
    for (int i = 0; i < nSynapses_; ++i)
        w_[i] = new_w[i];
    w_[ nSynapses_ ] = 0.0;
    w_[ nSynapses_ + 1] = -0.01;

    if (!allow_initial_inhibatory_weights)
        for (int i = 0; i < nSynapses_; ++i)
            if (w_[i] < 0.0)
                throw runtime_error(string("\nError Trying to set weights of Tempotron. Negative weight encountered while allow_inhibatory_weights=false\n"));
    if (restrict_weights_sign)
        enable_weight_sign_restriction();

}

void Tempotron2::rand_w(double std) {
    NR::GausianDist rnd(0.0, std);
    for (int i = 0; i < nSynapses_; w_[i++] = rnd());
    w_[ nSynapses_ ] = 0.0;
    w_[ nSynapses_ + 1] = -0.01;
    delta_w_gag_psik = vector<double>(nSynapses_ + 1, 0.0);

    if (!allow_initial_inhibatory_weights)
        for (int i = 0; i < nSynapses_; ++i)
            if (w_[i] < 0.0)
                w_[i] *= -1.0;
    if (restrict_weights_sign)
        enable_weight_sign_restriction();

#ifdef USE_PYTHON
     pyw_begin=w_.begin();
      pyw_end=w_.end()-2; 
#endif
}

void Tempotron2::traceVoltage(double dt) {
    static double V_next, d, d_s;
    deque<Spike2> grided = *spikes_;
    std::deque< Spike2 >::iterator spikesEnd, cSpike, nSpike;
    voltageTrace_ = deque< MaxPoint > ();

    for (double t = 0.0; t <= K_->T+3*K_->tau; t += dt)
        grided.push_back(Spike2(nSynapses_, t, exp(t / K_->tau), exp(t / K_->tau_s)));
    sort(grided.begin(), grided.end());

    spikesEnd = grided.end();
    cSpike = grided.begin();
    nSpike = ++(grided.begin());

    d = 0.0;
    d_s = 0.0;

    voltageTrace_.push_back(MaxPoint(0.0,0.0));
    while (nSpike != spikesEnd) {
        d += w_[ cSpike->affarent ] * cSpike->tau_exponent;
        d_s += w_[ cSpike->affarent ] * cSpike->tau_s_exponent;

        if (nSpike->affarent == nSynapses_) {
            V_next = K_->Vo * (d / nSpike->tau_exponent - d_s / nSpike->tau_s_exponent);
            voltageTrace_.push_back(MaxPoint(nSpike->time, V_next));
        }

        ++cSpike;
        ++nSpike;
    }
}

void Tempotron2::loadVector(string vectorFileName, ostream& out /*= cout*/) {
    ifstream vin(vectorFileName.c_str(), ios::binary);
    if (vin.is_open()) {
        char * data = new char[ sizeof (double) ];

        for (int i = 0; i < w_.size(); ++i) {
            vin.read(data, sizeof (double));
            w_[i] = *(reinterpret_cast<double *> (data));
        }

        vin.close();
        delta_w_gag_psik = vector<double>(nSynapses_, 0.0);
        out << "\nLoaded conection vector from file: " << vectorFileName << "\n";
    } else out << "ERROR: Tempotron2::loadVector() - file not open!" << std::endl;
}

void Tempotron2::saveVector(std::string vectorFN, std::ostream& out /*= std::cout*/) {
    std::ofstream vout(vectorFN.c_str(), ios::binary);
    if (vout.is_open()) {
        for (int i = 0; i < w_.size(); ++i)
            vout.write(reinterpret_cast<char *> (&(w_[i])), sizeof (double));
        vout.close();
        out << "\nSaved connection vector to : " << vectorFN << std::endl;
    } else out << "ERROR: Tempotron2::loadVector() - file not open!" << std::endl;
}

void Tempotron2::serializeVectors(std::ostream& o) {
    int s = w_.size();
    o.write(reinterpret_cast<char *> (&s), sizeof (int));
    for (int i = 0; i < w_.size(); ++i)
        o.write(reinterpret_cast<char *> (&(w_[i])), sizeof (double));
    s = delta_w_gag_psik.size();
    o.write(reinterpret_cast<char *> (&s), sizeof (int));
    for (int i = 0; i < delta_w_gag_psik.size(); ++i)
        o.write(reinterpret_cast<char *> (&(delta_w_gag_psik[i])), sizeof (double));
    if (restrict_weights_sign) {
        s = w_sign_.size();
        o.write(reinterpret_cast<char *> (&s), sizeof (int));
        for (int i = 0; i < w_sign_.size(); ++i)
            o.write(reinterpret_cast<char *> (&(w_sign_[i])), sizeof (double));
    }

}

void Tempotron2::loadVectors(std::istream& i) {
    char int_data[sizeof (int) ];
    char double_data[sizeof (double) ];

    int s;
    cout << "loading weights" << endl;
    i.read(int_data, sizeof (int));
    s = *(reinterpret_cast<int *> (int_data));
    w_ = vector<double>(s, 0.0);

    for (int j = 0; j < w_.size(); ++j) {
        i.read(double_data, sizeof (double));
        w_[j] = *(reinterpret_cast<double *> (double_data));
    }

    i.read(int_data, sizeof (int));
    s = *(reinterpret_cast<int *> (int_data));
    delta_w_gag_psik = vector<double>(s, 0.0);
    cout << "loaing delta_w_gag_psik" << endl;
    for (int j = 0; j < delta_w_gag_psik.size(); ++j) {
        i.read(double_data, sizeof (double));
        delta_w_gag_psik[j] = *(reinterpret_cast<double *> (double_data));
    }
    if (restrict_weights_sign) {
        i.read(int_data, sizeof (int));
        s = *(reinterpret_cast<int *> (int_data));
        w_sign_ = vector<double>(s, 0.0);

        for (int j = 0; j < w_sign_.size(); ++j) {
            i.read(double_data, sizeof (double));
            w_sign_[j] = *(reinterpret_cast<double *> (double_data));
        }
    }



}

double Tempotron2::KminusTh::operator()(long double t) const {
    return ( tmp->K_->Vo * (tmp->D * std::exp(-t / tmp->K_->tau) - tmp->D_s * std::exp(-t / tmp->K_->tau_s)) / tmp->threshold_-1.);
}

void Tempotron2::weightDecay(double lambda) {
    for (int i=0; i<max_aff_learn; ++i)
        w_[i] *= lambda;
}

