/* 
 * File:   RanUtilsRandom.h
 * Author: ranr
 *
 * Created on March 30, 2009, 10:37 PM
 */
// THIS IS A PARTIAL HEADER FILE COMPARED TO THE LIBRANUTILS PROJECT 
#ifndef _RANUTILSRANDOM_H
#define	_RANUTILSRANDOM_H

#include "NR.h"
#include <vector>

namespace RanUtils {

    class Random {
    private:
        NR::Ran * rnd;
    public:

        Random() {
            rnd = &(NR::Rnd);
        }

        Random(NR::Ran & ranGen) {
            rnd = &ranGen;
        }

        inline int operator() () {
            return rnd->int32();
        }

        int operator() (int max) {
            return (int) (max * rnd->doub());
        }

        int nextInt(int max) {// returns uniform random integer in [0,max)
            return (int) (max * rnd->doub());
        }

        int nextInt(int min, int max) {// returns uniform random integer in [min,max)
            return min + nextInt(max - min);
        }

        float nextFloat() {
            return (float) (rnd->doub());
        }

        double nextDouble() {
            return rnd->doub();
        }

        float nextFloat(float max) {
            return max * (float) (rnd->doub());
        }

        float nextFloat(float min, float max) {
            return min + (max - min)*(rnd->doub());
        }

        double nextDouble(double max) {
            return max * (rnd->doub());
        }

        double nextDouble(double min, double max) {
            return min + (max - min)*(rnd->doub());
        }

        bool nextBool(double f) {
            return ( rnd->doub() < f) ? true : false;
        }

        bool nextBool() {
            return nextBool(0.5);
        }

        void printID(std::ostream& o = std::cout) {
            o << " rnd add.: " << rnd << " seed = " << (*rnd).getSeed() << std::endl;
        }
    };

    
    class PoissonProccess {//this is a poisson process dist...
    public:
        
        PoissonProccess() {
            tau_ = 1.0;
            rnd_ = &(NR::rnd2);
            refractoryPeriod=0.0;
        }

        PoissonProccess(double tau) {
            tau_ = tau;
            rnd_ = &(NR::rnd2);
            refractoryPeriod=0.0;
        }

        inline double NextEvent() {
            static double x;
            do
                x = (*rnd_)(); while (x == 0);
            return refractoryPeriod-std::log(x) * tau_;
        }

        std::vector<double> NextEventsInTimeT(double T) {
            std::vector<double> ans;
            double tt = 0.0;

            while (tt < T) {
                tt += this->NextEvent();
                ans.push_back(tt);
            }
            ans.pop_back();

            return ans;
        }

        std::vector<double> NextNEvents(int N) {
            std::vector<double> ans(N);
            double t = 0;

            for (int i = 0; i < N; ++i) {
                t += this->NextEvent();
                ans[i] = t;
            }
            return ans;
        }

        void setTau(double tau) {
            tau_ = tau;
        }

        double getTau() {
            return tau_;
        }

        void setRanGen(NR::Ran * rnd) {
            rnd_ = rnd;
        }

        void SetRefractoryPeriod(double refractoryPeriod) {
            this->refractoryPeriod = refractoryPeriod;
        }

        double GetRefractoryPeriod() const {
            return refractoryPeriod;
        }
    protected:
        double tau_,refractoryPeriod;
        NR::Ran* rnd_;

    };

 typedef PoissonProccess PoissonDistribution;
}

namespace RU = RanUtils;



#endif	/* _RANUTILSRANDOM_H */

