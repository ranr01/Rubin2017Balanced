#ifndef NRC_H
#define NRC_H

#include <ctime>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <string>
#include <sstream>
#ifdef CPP11
        #define CONSTEXPR constexpr
#else
        #define CONSTEXPR const
#endif
namespace NR {

    void nrerror(const char * error_text);
    
    class ran2
    {
    protected:
        static const long IM1 = 2147483563;
        static const long IM2 = 2147483399;
        static CONSTEXPR double AM = (1.0/IM1);
        static const long IMM1 = (IM1-1);
        static const int IA1 = 40014;
        static const int IA2 = 40692;
        static const int IQ1 = 53668;
        static const int IQ2 = 52774;
        static const int IR1 = 12211;
        static const int IR2 = 3791;
        static const int NTAB = 32;
        static const long NDIV = (1+IMM1/NTAB);
        static CONSTEXPR double EPS = 1.2e-7;
        static CONSTEXPR double RNMX = 0.99999988;
        long idum;
        long seed_;


    public:
        ran2(long seed) { seed_=0; setSeed(seed);}
        ran2() { seed_ =  -time(NULL); idum = seed_;}

        long getSeed() { return seed_;}
        void setSeed(long seed) { 
            if (seed<0){
                seed_ = seed; 
                idum = seed_;
                std::cout<<"\nSeed set for ran2. seed="<<seed_<<std::endl;
            } else
              std::cout<<"\nERROR: seed must be negative."<<std::endl;
            if (seed_==0){
                seed_ =  -time(NULL); 
                idum = seed_;
            }
        }

        double operator()();
    };

    class Ran {
/*Implementation of the highest quality recommended generator. The constructor is called with
an integer seed and creates an instance of the generator. The member functions int64, doub,
and int32 return the next values in the random sequence, */
    protected:
        unsigned long u,v,w;
        unsigned long seed_;
    public:
        Ran(unsigned long j) : screen_output(true) {
        //Constructor. Call with any integer seed (except value of v in setSeed).
          setSeed(j);
        }
        Ran() : screen_output(true)  { setSeed((unsigned long)(time(NULL))); }

        //Backword compatibility constructor for negative seed
        Ran(long j) : screen_output(true)  {
            if (j>0) setSeed((unsigned long)(j));
            else setSeed((unsigned long)(-j));
        }
        unsigned long getSeed() { return seed_;}
        void setSeed(unsigned long j){
           seed_=j;
           v=4101842887655102017LL; w=1;
           if (j!=v){
           u = j ^ v; int64();
           v = u; int64();
           w = v; int64();
           if (screen_output)
               std::cout<<"Ran::setSeed : RNG address - "<<this<<" seed has been set. seed_="<<seed_<<std::endl;
           } else {
               std::cout<<"Ran::setSeed : ERROR ileagal seed. seting seed to time(NULL)\n";
               setSeed((unsigned long)(time(NULL)));
           }
        }
        void setSeed(long j){
           if (j>0) setSeed((unsigned long)(j));
            else setSeed((unsigned long)(-j));
        }
        void setSeed(int j){
           if (j>0) setSeed((unsigned long)(j));
            else setSeed((unsigned long)(-j));
        }

        inline unsigned long int64() { //  Return 64-bit random integer. See text for explanation of method.
            u = u * 2862933555777941757LL + 7046029254386353087LL;
            v ^= v >> 17; v ^= v << 31; v ^= v >> 8;
            w = 4294957665U*(w & 0xffffffff) + (w >> 32);
            unsigned long x = u ^ (u << 21); x ^= x >> 35; x ^= x << 4;
            return (x + v) ^ w;
        }

        inline double doub() { return 5.42101086242752217E-20 * int64(); }
        inline double operator()() {return 5.42101086242752217E-20 * int64();}
        //Return random double-precision floating value in the range 0. to 1.
        inline int int32() { return (int)(int64()); }
        // Return 32-bit random integer.

        bool screen_output;
    };


    extern Ran Rnd;
    //Backward comptibility to programs using ran2. The code works the same way but uses Ran instead.
    extern Ran & rnd2;

    class GausianDist
    {
    private:
        float mean_, std_;
        double gset,fac,rsq,v1,v2;
        bool iset;
        Ran* rnd;

    public:
        GausianDist(float mean, float std, Ran & ranGen ) : mean_(mean), std_(std), iset(false) { rnd = &ranGen;}
        GausianDist(float mean, float std) : mean_(mean), std_(std), iset(false) { rnd = &Rnd;}
        GausianDist() : mean_(0.0), std_(1.0) , iset(false) { rnd = &Rnd;}
        double getMean() {return mean_;}
        double getStd() {return std_;}
        
        double operator()();
    };


    namespace Math
    {
        float gammln(float xx); //Returns the value ln[ (xx)] for xx > 0. 
        float gammp(float a, float x);
        float gammq(float a, float x);
        void gcf(float *gammcf, float a, float x, float *gln);
        void gser(float *gamser, float a, float x, float *gln);
        float erf(float x);
        float erfc(float x);
        double SIGN(double a,double b);
	void gauher(double x[], double w[], int n); //returns weigts and abscissas for gaussion quadrature woth gauss hermite integration

        template<class func> 
        double find_x(double x1, double x2, func const & f, double xacc = 1e-10,double facc = 1e-8)  throw(std::runtime_error)
        {
            const int MAXIT =  6000;
            const double UNUSED = (-1.11e30);

            int j;
            double ans,fh,fl,fm,fnew,s,xh,xl,xm,xnew;

            fl=f(x1);
            fh=f(x2);
            //std::cout<<"f("<<x1<<")="<<fl<<" f("<<x2<<")="<<fh<<" "<<std::flush;
            if ((fl > 0.0 && fh < 0.0) || (fl < 0.0 && fh > 0.0))
            {
                xl=x1;
                xh=x2;
                ans=UNUSED;
                for (j=1;j<=MAXIT;j++)
                {
                    xm=0.5*(xl+xh);
                    fm=f(xm);
                    s=std::sqrt(fm*fm-fl*fh);
                    if (s == 0.0)
                        return ans;
                    xnew=xm+(xm-xl)*(( (fl >= fh) ? 1.0 : -1.0 )*fm/s);
                    if ( (std::fabs(xnew-ans) <= xacc) && (std::fabs(fm) <= facc) )
                        return ans;
                    ans=xnew;
                    fnew=f(ans);
                    if (fnew == 0.0)
                        return ans;
                    if (SIGN(fm,fnew) != fm)
                    {
                        xl=xm;
                        fl=fm;
                        xh=ans;
                        fh=fnew;
                    }
                    else if (SIGN(fl,fnew) != fl)
                    {
                        xh=ans;
                        fh=fnew;
                    }
                    else if (SIGN(fh,fnew) != fh)
                    {
                        xl=ans;
                        fl=fnew;
                    }
                    else throw std::runtime_error(std::string("NR::MATH::find_x ERROR : never get here."));
                    if ( (std::fabs(xh-xl) <= xacc) && (std::fabs(fnew) <= facc) )
                        return ans;
                }
                std::cout<<"xh-xl="<<xh-xl<<" f("<<xl<<")="<<fl<<" f("<<xh<<")="<<fh<<" "<<std::flush;
                throw std::runtime_error(std::string("NR::MATH::find_x ERROR : zriddr exceed maximum iterations."));
            }
            else
            {
                if (fl == 0.0) return x1;
                if (fh == 0.0) return x2;
                std::ostringstream o;
                o<<"NR::MATH::find_x ERROR : root must be bracketed in zriddr. x1="<<x1<<" x2="<<x2<<" f(x1)="<<fl<<" f(x2)="<<fh;
                throw std::runtime_error(o.str());
            }
            return 0.0;
        }
    }



}

#endif
