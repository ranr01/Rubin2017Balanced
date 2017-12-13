#include <math.h>
#include <ctime>
#include "NR.h"
#include "NRfitting.h"


using namespace std;

double NR::chi_sqr(vector<double>& x,vector<double>& y,vector<double>& grad,vector< vector<double> >& hess,NR::model_function & f)
{
  double chi=0.0;
  double tmpf,del;
  vector<double> tmpgrad;
  grad=vector<double>(f.ma,0.0);
  hess=vector< vector<double> >(f.ma,vector<double>(f.ma,0.0));

  for (int i=0 ; i<x.size(); ++i)
    {
      tmpf=f(x[i]);
      del=(y[i]-tmpf);
      chi+=del*del;
      tmpgrad=f.grad(x[i]);
      for(int k=0 ; k<f.ma ; ++k){
	grad[k]-=2*del*tmpgrad[k];
	for(int l=k; l<f.ma ; ++l) hess[k][l]+=2*tmpgrad[k]*tmpgrad[l];
      }
    }

  for(int k=0; k<f.ma ; ++k)
    for(int l=k ; l<f.ma ; ++l)
      hess[l][k]=hess[k][l];

  return chi;

}

void NR::mrqmin2D(vector<double>& x,vector<double>& y, NR::model_function& f)
{
  bool stop=false;
  vector<double> grad, old_grad,old_a;
  vector< vector<double> > hess, old_hess;
  static double max_iter=10000;
  double current_chi_sqr,old_chi_sqr,lambda;

  lambda=0.001;
  current_chi_sqr=chi_sqr(x,y,grad,hess,f);

  int count=0;

  while ((!stop)&&(count<max_iter)){
    /* 2d solution */
    double det=hess[0][0]*hess[1][1]*(1+lambda)*(1+lambda)-hess[0][1]*hess[0][1];
    double da0=-(hess[1][1]*(1+lambda)*grad[0]-hess[0][1]*grad[1])/det;
    double da1=-(hess[0][0]*(1+lambda)*grad[1]-hess[0][1]*grad[0])/det;
    /***************/
    old_a=f.a;
    f.a[0]+=da0;
    f.a[1]+=da1;
    old_chi_sqr=current_chi_sqr;
    old_grad=grad;
    old_hess=hess;

    current_chi_sqr=chi_sqr(x,y,grad,hess,f);
    
    if (current_chi_sqr < old_chi_sqr) {
      stop= ( (1.0-current_chi_sqr/old_chi_sqr)<1e-4 );
      lambda *= 0.1;
    } else {
      lambda *= 10.0;
      current_chi_sqr=old_chi_sqr;
      f.a=old_a;
      grad=old_grad;
      hess=old_hess;
    }
    ++count;
  }
}


void NR::nrerror(const char * error_text){
    std::cout<<error_text<<std::endl;
}

float NR::Math::gammln(float xx) //Returns the value ln[ (xx)] for xx > 0. 
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}

float NR::Math::gammp(float a, float x)
{
    float gamser,gammcf,gln;

	if (x < 0.0 || a <= 0.0) NR::nrerror("Invalid arguments in routine gammp");
	if (x < (a+1.0)) {
		NR::Math::gser(&gamser,a,x,&gln);
		return gamser;
	} else {
		NR::Math::gcf(&gammcf,a,x,&gln);
		return 1.0-gammcf;
	}
}

float NR::Math::gammq(float a, float x)
{
    float gamser,gammcf,gln;

	if (x < 0.0 || a <= 0.0) NR::nrerror("Invalid arguments in routine gammq");
	if (x < (a+1.0)) {
		NR::Math::gser(&gamser,a,x,&gln);
		return 1.0-gamser;
	} else {
		NR::Math::gcf(&gammcf,a,x,&gln);
		return gammcf;
	}
}

void NR::Math::gcf(float *gammcf, float a, float x, float *gln)
{
	const int itmax = 100;
    const float eps = 3.0e-7;
    const float fpmin = 1.0e-30;
	int i;
	float an,b,c,d,del,h;

	*gln=NR::Math::gammln(a);
	b=x+1.0-a;
	c=1.0/fpmin;
	d=1.0/b;
	h=d;
	for (i=1;i<=itmax;i++) {
		an = -i*(i-a);
		b += 2.0;
		d=an*d+b;
		if (fabs(d) < fpmin) d=fpmin;
		c=b+an/c;
		if (fabs(c) < fpmin) c=fpmin;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) < eps) break;
	}
	if (i > itmax) NR::nrerror("a too large, itmax too small in gcf");
	*gammcf=exp(-x+a*log(x)-(*gln))*h;
}

void NR::Math::gser(float *gamser, float a, float x, float *gln)
{
    const int itmax = 100;
    const float eps = 3.0e-7;
    int n;
	float sum,del,ap;

	*gln=NR::Math::gammln(a);
	if (x <= 0.0) {
		if (x < 0.0) NR::nrerror("x less than 0 in routine gser");
		*gamser=0.0;
    	return;
	} else {
		ap=a;
		del=sum=1.0/a;
		for (n=1;n<=itmax;n++) {
			++ap;
			del *= x/ap;
			sum += del;
			if (fabs(del) < fabs(sum)*eps) {
				*gamser=sum*exp(-x+a*log(x)-(*gln));
				return;
			}
		}
		NR::nrerror("a too large, itmax too small in routine gser");
		return;
	}
}

float NR::Math::erf(float x)
{
    return x < 0.0 ? -NR::Math::gammp(0.5,x*x) : NR::Math::gammp(0.5,x*x);
}

float NR::Math::erfc(float x)
{
    return x < 0.0 ? 1.0+NR::Math::gammp(0.5,x*x) : NR::Math::gammq(0.5,x*x);
}


 NR::Ran NR::Rnd;
 NR::Ran & NR::rnd2(NR::Rnd);

 double NR::Math::SIGN(double a,double b) { return ( (b >= 0.0)? std::fabs(a) : -std::fabs(a) ); }

 double NR::ran2::operator()()
   { 
       int j;
       long k;
       static long idum2=123456789;
       static long iy=0;
       static long iv[NTAB];
       float temp;

       if (idum <= 0) {
           if (-(idum) < 1) idum=1;
           else idum = -(idum);
           idum2=(idum);
           for (j=NTAB+7;j>=0;--j) {
               k=(idum)/IQ1;
               idum=IA1*(idum-k*IQ1)-k*IR1;
               if (idum < 0) idum += IM1;
               if (j < NTAB) iv[j] = idum;
           }
           iy=iv[0];
       }
       k=(idum)/IQ1;
       idum=IA1*(idum-k*IQ1)-k*IR1;
       if (idum < 0) idum += IM1;
       k=idum2/IQ2;
       idum2=IA2*(idum2-k*IQ2)-k*IR2;
       if (idum2 < 0) idum2 += IM2;
       j=iy/NDIV;
       iy=iv[j]-idum2;
       iv[j] = idum;
       if (iy < 1) iy += IMM1;
       if ((temp=AM*iy) > RNMX) return RNMX;
       else return temp;
 }

 double NR::GausianDist::operator() () 
      {

         static double u,v,x,y,q;
         do {
            u = (*rnd)();
            v = 1.7156*((*rnd)()-0.5);
            x = u - 0.449871;
            y = std::abs(v) + 0.386595;
            q = x*x+ y*(0.19600*y-0.25472*x);
         } while (q > 0.27597
                       && (q > 0.27846 || v*v > -4.*std::log(u)*u*u));
        return mean_ + std_*v/u;

        /*  if (!iset)
          { //We don t have an extra deviate handy, so 
              do { 
                  v1=2.0*(*rnd)()-1.0; //pick two uniform numbers in the square extending from -1 to +1 in each direction, 
                  v2=2.0*(*rnd)()-1.0; 
                  rsq=v1*v1+v2*v2;   //see if they are in the unit circle,and if they are not, try again. 
              } while (rsq >= 1.0 || rsq == 0.0);  
              
              //Now make the Box-Muller transformation to get two normal deviates. 
              //Return one and save the other for next time. 
              fac=sqrt(-2.0*std::log(rsq)/rsq); 
              gset=v1*fac; 
              iset=true; //Set flag. 
              return (std_*v2*fac + mean_); 
          } 
          else 
          { //We have an extra deviate handy,so unset the flag, and return it
              iset=false;  
              return (std_*gset + mean_); 
          }*/
        } 

void NR::Math::gauher(double x[], double w[], int n)
{
        static double EPS=3.0e-14 , PIM4=0.7511255444649425;
	static int MAXIT=10;

	int i,its,j,m;
	double p1,p2,p3,pp,z,z1;

	m=(n+1)/2;
	for (i=1;i<=m;i++) {
		if (i == 1) {
			z=sqrt((double)(2*n+1))-1.85575*pow((double)(2*n+1),-0.16667);
		} else if (i == 2) {
			z -= 1.14*pow((double)n,0.426)/z;
		} else if (i == 3) {
			z=1.86*z-0.86*x[1];
		} else if (i == 4) {
			z=1.91*z-0.91*x[2];
		} else {
			z=2.0*z-x[i-2];
		}
		for (its=1;its<=MAXIT;its++) {
			p1=PIM4;
			p2=0.0;
			for (j=1;j<=n;j++) {
				p3=p2;
				p2=p1;
				p1=z*sqrt(2.0/j)*p2-sqrt(((double)(j-1))/j)*p3;
			}
			pp=sqrt((double)2*n)*p2;
			z1=z;
			z=z1-p1/pp;
			if (fabs(z-z1) <= EPS) break;
		}
		if (its > MAXIT) nrerror("too many iterations in gauher");
		x[i]=z;
		x[n+1-i] = -z;
		w[i]=2.0/(pp*pp);
		w[n+1-i]=w[i];
	}
}
