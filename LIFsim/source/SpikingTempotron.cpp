//#include <iostream>
#include <deque>

#include <numeric>
#include <stdexcept>
#include "SpikingTempotron.h"
#include "Spike2.h"
#include "Tempotron2.h"
#include "NR.h"
//#include "InputLayer2.h"

using namespace std;

void SpikingTempotron::activate() {
    //removed from this code since no learning is needed
}

void SpikingTempotron::activate_noTeacher() {
    static long double lnD, lnD_s, t_max, V_cur, V_next;
    static std::deque< Spike2 >::const_iterator spikesEnd;
    static int dtb;

    KminusTh f(this);
    crossings_.clear();

#ifdef DEBUG
    int stage;
#endif
    spikesEnd = spikes_->end();

    currentSpike_ = spikes_->begin();
    nextSpike_ = currentSpike_;
    ++nextSpike_;

    timeBlock_ = 0;
    expTimeBlock_taum_ = std::exp((long double) (-spikes_->timeBlockSize / K_->tau));
    expTimeBlock_taus_ = std::exp((long double) (-spikes_->timeBlockSize / K_->tau_s));
    TB_BeginTime_ = 0.0;

    reset_D_ = 0.0;
    D = 0.0;
    D_s = 0.0;

    try {
        while (nextSpike_ != spikesEnd) {
            if (currentSpike_->timeBlock > timeBlock_) {
                D *= std::pow(expTimeBlock_taum_, currentSpike_->timeBlock - timeBlock_);
                reset_D_ *= std::pow(expTimeBlock_taum_, currentSpike_->timeBlock - timeBlock_);
                D_s *= std::pow(expTimeBlock_taus_, currentSpike_->timeBlock - timeBlock_);
                timeBlock_ = currentSpike_->timeBlock;
                TB_BeginTime_ = timeBlock_ * spikes_->timeBlockSize;
            }

            D += w_[ currentSpike_->affarent ] * currentSpike_->tau_exponent;
            D_s += w_[ currentSpike_->affarent ] * currentSpike_->tau_s_exponent;

            if ((D > 0.0) && (D_s > 0.0)) {
                lnD = std::log(D);
                lnD_s = std::log(D_s);
                t_max = K_->t_max_fac * (lnD_s - lnD + K_->lntt);
                if (t_max > currentSpike_->time) {
                    dtb = nextSpike_->timeBlock - currentSpike_->timeBlock;
                    if (t_max <= nextSpike_->time + dtb * spikes_->timeBlockSize) {
                        V_max_ = std::exp(K_->v_max_fac * (K_->tau * lnD - K_->tau_s * lnD_s));
                        if (V_max_ > threshold_) {
#ifdef DEBUG
                            stage = 1;
#endif
                            crossings_.push_back(TB_BeginTime_ + NR::Math::find_x(currentSpike_->time, t_max, f));
                            reset_V(crossings_.back());
                        }
                    } else {
                        V_next = K_->Vo * (D / nextSpike_->tau_exponent * std::pow(expTimeBlock_taum_, dtb) -
                                D_s / nextSpike_->tau_s_exponent * std::pow(expTimeBlock_taus_, dtb));
                        if (V_next > threshold_) {
#ifdef DEBUG
                            stage = 2;
#endif
                            crossings_.push_back(TB_BeginTime_ + NR::Math::find_x(currentSpike_->time, nextSpike_->time + dtb * spikes_->timeBlockSize, f));
                            reset_V(crossings_.back());
                        }
                    }
                }
            }
            ++currentSpike_;
            ++nextSpike_;
        }
        if (currentSpike_->timeBlock > timeBlock_) {
            D *= std::pow(expTimeBlock_taum_, currentSpike_->timeBlock - timeBlock_);
            reset_D_ *= std::pow(expTimeBlock_taum_, currentSpike_->timeBlock - timeBlock_);
            D_s *= std::pow(expTimeBlock_taus_, currentSpike_->timeBlock - timeBlock_);
            timeBlock_ = currentSpike_->timeBlock;
            TB_BeginTime_ = timeBlock_ * spikes_->timeBlockSize;
        }

        D += w_[ currentSpike_->affarent ] * currentSpike_->tau_exponent;
        D_s += w_[ currentSpike_->affarent ] * currentSpike_->tau_s_exponent;

        if ((D > 0.0) && (D_s > 0.0)) {
            lnD = std::log(D);
            lnD_s = std::log(D_s);
            t_max = K_->t_max_fac * (lnD_s - lnD + K_->lntt);

            if (t_max > currentSpike_->time) {
                V_max_ = std::exp(K_->v_max_fac * (K_->tau * lnD - K_->tau_s * lnD_s));
                if (V_max_ > threshold_) {
#ifdef DEBUG
                    stage = 7;
#endif
                    crossings_.push_back(TB_BeginTime_ + NR::Math::find_x(currentSpike_->time, t_max, f));
                    reset_V_last_spike(crossings_.back());
                }
            }
        }

    } catch (std::runtime_error & e) {
        cout << " Error in activate_noTeacher()\n";
#ifdef DEBUG
        cout << "stage=" << stage;
#endif
        cout << " S" << currentSpike_ - spikes_->begin() << " " << D << " " << D_s << "\n"
                << e.what() << "\n"
                << " t_max=" << t_max << " V_next=" << V_next << " V_max_=" << V_max_ << " thershold_=" << threshold_ << "\n";
        cout << "curruntSpike: " << *currentSpike_ << '\n'
                << "nextSpike: " << *nextSpike_ << "\n"
                << "End=" << (nextSpike_ == spikesEnd) << endl;
        throw;
    }

    ++currentSpike_;

  
}

int SpikingTempotron::learn() {
    //Removed from this code since no learning is needed  
    return 0;
}

bool SpikingTempotron::td_err(int count) {
    static long double V;
    if (count) {
        V = (D_td_[count] / spikes_->tds[count]->tau_exponent - D_s_td_[count] / spikes_->tds[count]->tau_s_exponent);
        V -= (D_td_[count - 1] / spikes_->tds[count - 1]->tau_exponent - D_s_td_[count - 1] / spikes_->tds[count - 1]->tau_s_exponent)
                * spikes_->tds[count - 1]->tau_exponent / spikes_->tds[count]->tau_exponent;

        return ((V * K_->Vo) < threshold_);
    } else
        return ((K_->Vo * (D_td_[count] / spikes_->tds[count]->tau_exponent - D_s_td_[count] / spikes_->tds[count]->tau_s_exponent)) < threshold_);

}

void SpikingTempotron::reset_V(double t) {
    static long double lnD, lnD_s, t_max, V_next, dD;
    static int dtb;
    KminusTh f(this);
#ifdef DEBUG
    static int spike_count = 0;
#endif

    t -= TB_BeginTime_;
    dD = (long double)(fraction_reset*threshold_ / K_->Vo) * std::exp((long double) (t / K_->tau));
    D -= dD;
    reset_D_ -= dD;

    if ((D > 0.0) && (D_s > 0.0)) {
        lnD = std::log(D);
        lnD_s = std::log(D_s);
        t_max = K_->t_max_fac * (lnD_s - lnD + K_->lntt);
        if (t_max > t) {
            dtb = nextSpike_->timeBlock - timeBlock_;
            if (t_max <= nextSpike_->time + dtb * spikes_->timeBlockSize) {
                V_max_ = std::exp(K_->v_max_fac * (K_->tau * lnD - K_->tau_s * lnD_s));
                if (V_max_ > threshold_) {
#ifdef DEBUG
                    ++spike_count;
                    long double V_t = K_->Vo * ( D * std::exp((long double)(-t/K_->tau)) - D_s * std::exp((long double)(-t/K_->tau_s)));
                    cout<<"Multiple Spikes (t_max): t="<<t<<" V(t)="<<V_t<<" t_max="<<t_max<<" V_max="<<V_max_<<" spike_count="<<spike_count<<endl;
#endif
                    crossings_.push_back(TB_BeginTime_ + NR::Math::find_x(t, t_max, f));
                    reset_V(crossings_.back());
                }
            } else {
                V_next = K_->Vo * (D / nextSpike_->tau_exponent * std::pow(expTimeBlock_taum_, dtb) -
                        D_s / nextSpike_->tau_s_exponent * std::pow(expTimeBlock_taus_, dtb));
                if (V_next > threshold_) {
#ifdef DEBUG
                    ++spike_count;
                    long double V_t = K_->Vo * ( D * std::exp((long double)(-t/K_->tau)) - D_s * std::exp((long double)(-t/K_->tau_s)));
                    cout<<"Multiple Spikes (t_next): t="<<t<<" V(t)="<<V_t<<" t_next="<<nextSpike_->time + dtb * spikes_->timeBlockSize
                            <<" V_next="<<V_next<<" spike_count="<<spike_count<<endl;
#endif

                    crossings_.push_back(TB_BeginTime_ + NR::Math::find_x(t, nextSpike_->time + dtb * spikes_->timeBlockSize, f));
                    reset_V(crossings_.back());
                }
            }
        }
    }
#ifdef DEBUG
    spike_count=0;
#endif
}
    void SpikingTempotron::reset_V_last_spike(double t) {
    static long double lnD, lnD_s, t_max, dD;
    KminusTh f(this);

    t -= TB_BeginTime_;
    dD = (fraction_reset*threshold_ / K_->Vo) * std::exp((long double) (t / K_->tau));
    D -= dD;
    reset_D_ -= dD;

    if ((D > 0.0) && (D_s > 0.0)) {
        lnD = std::log(D);
        lnD_s = std::log(D_s);
        t_max = K_->t_max_fac * (lnD_s - lnD + K_->lntt);
        if (t_max > t) {
                V_max_ = std::exp(K_->v_max_fac * (K_->tau * lnD - K_->tau_s * lnD_s));
                if (V_max_ > threshold_) {
                    crossings_.push_back(TB_BeginTime_ + NR::Math::find_x(t, t_max, f));
                    reset_V_last_spike(crossings_.back());
                }
            }
    }
}

void SpikingTempotron::teacher_reset_V(long double Vtd) {

    static long double dD;
    D_td_[td_count_] = D - reset_D_;
    D_s_td_[td_count_] = D_s;
    V_td_[td_count_] = Vtd;

    dD = (Vtd / K_->Vo) * this->tdSpike_->tau_exponent;
    D -= dD;
    reset_D_ -= dD;
}

