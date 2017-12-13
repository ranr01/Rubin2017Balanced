/* 
 * File:   SpikeTrain.h
 * Author: ranr
 *
 * Created on June 11, 2013, 10:04 PM
 */

#ifndef SPIKETRAIN_H
#define	SPIKETRAIN_H
#include  "Spike2.h"
#include <deque>
#include <iostream>
#include "PSP_Kernel2.h"
#include <vector>
#include <algorithm>
#include <stdexcept>

class SpikeTrain : public std::deque< Spike2 > {
public:
    SpikeTrain() : std::deque< Spike2 >(), number_of_tds(0) , timeBlockSize(1.0){}
    virtual ~SpikeTrain() { }

    SpikeTrain(SpikeTrain::iterator b,SpikeTrain::iterator e) : std::deque< Spike2 >(b,e), number_of_tds(0)  , 
            timeBlockSize(1.0){}
    
    void clear() {
        //std::cout<<"*"<<std::flush;
        number_of_tds=0;
        tds.clear();
        tds_minus_eps.clear();
        ((std::deque< Spike2 > *)(this))->clear();
    }
    void add_spike(int aff, double time){
        push_back(Spike2(aff,time));
    }
    void sort(){
        std::sort(begin(),end());
    }
    //void toFile(std::ostream & o);

    void applyTimeBlock(){
        int tb;
        for (SpikeTrain::iterator it=begin(); it!=end(); ++it){
            tb = (int)(it->time / timeBlockSize);
            it->time -= timeBlockSize*tb;
            it->timeBlock = tb;
        }
    }
    void calcExponents(PSP_Kernel2 const & K){
        for (SpikeTrain::iterator it=begin(); it!=end(); ++it){
            it->tau_exponent = std::exp((long double)(it->time/K.tau));
            it->tau_s_exponent = std::exp((long double)(it->time/K.tau_s));
        }
    }
    
    void setTdMarkers(int N,bool tdOnly = false){
        // creating pointers to markers
        if (tdOnly){
                 tds = std::vector< std::deque<Spike2>::const_iterator > ();
                 tds_minus_eps = std::vector< std::deque<Spike2>::const_iterator > ();
     
                for (SpikeTrain::const_iterator it = begin(); (it != end()); ++it) {
                        if (it->affarent > N){
                                throw std::runtime_error(std::string("N != number of afferents"));
                        }
                        if (it->affarent == N) 
                            tds.push_back(it);
                }
                number_of_tds = tds.size();
        } else {
                bool inbetween = false;
        tds = std::vector< std::deque<Spike2>::const_iterator > ();
        tds_minus_eps = std::vector< std::deque<Spike2>::const_iterator > ();

      
          for (SpikeTrain::const_iterator it = begin(); (it != end()); ++it) {
              if (it->affarent > N){
                   throw std::runtime_error(std::string("N != number of afferents"));
              }
             if (it->affarent == N) {
                if (inbetween) {
                    tds.push_back(it);
                } else{
                    tds_minus_eps.push_back(it);
                }
                inbetween = !(inbetween);
            }
        }
        number_of_tds = tds.size();
        if  (tds.size()!=tds_minus_eps.size()){
            throw std::runtime_error(std::string("Non equal number of td and td-eps markers"));
        }
        }
    }
    
    int number_of_tds;
    double timeBlockSize;
    std::vector< SpikeTrain::const_iterator > tds, tds_minus_eps;
#ifdef USE_PYTHON
    std::vector<double> * getTds();
    std::vector<double> * getTdsMinusEps();
    
    bool operator==(SpikeTrain const & other){ return false; }
#endif
};


#endif	/* SPIKETRAIN_H */

