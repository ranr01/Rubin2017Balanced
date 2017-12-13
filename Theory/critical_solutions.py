import scipy.special
import scipy.optimize
import numpy as np
from numpy import *

# Eq. S48
H = lambda x: 0.5*scipy.special.erfc(x/np.sqrt(2.))
# Derivative of H
H_prime = lambda x: -exp(-x**2/2.)/sqrt(2*pi)

# Eq. S84
gamma= lambda x,s: 0.5*((x*x+1.)*H(s*x)+s*x*H_prime(x))
# Eq. S86
gamma_prime = lambda x,s: x*H(s*x)+s*H_prime(x)

# Solution to Eq. S94
@vectorize
def Delta(k,p_out):
    f = lambda D: p_out*(H_prime(D+k)-(D+k)*H(-D-k))-\
                (1.-p_out)*(H_prime(D-k)+(D-k)*H(D-k))
    return scipy.optimize.brentq(f,-10,10)

# Aux. function
J = lambda x: (x**2+1.)*H(x)+x*H_prime(x)
# Eq. S90
alpha_unconst = lambda D,k,p_out: 1./(p_out*J(-D-k)+(1-p_out)*J(D-k))

# Eq. S120 for the case of kappa_in
class Z_kappa_in(object):
    def __init__(self,kappa_in,p_out):
        self.kappa_in=kappa_in
        self.p_out=p_out

    def __call__(self,Q,llambda):
        return 1./(1.+(1.-Q/llambda**2)*self.G_Q(Q))

    # Eq. S121
    def G_Q(self,Q):
        k=self.kappa_in/sqrt(Q)
        D=Delta(k,self.p_out)
        return k*alpha_unconst(D,k,self.p_out)*2.*(1.-self.p_out)*\
              (H_prime(D-k)+(D-k)*H(D-k))


# R.H.S of eq. S92
def Q_eq(f_exc,B,phi,llambda,Q,Z):
    a = f_exc*gamma(B,1.)*Z(Q,1.)**2
    b = (1.-f_exc)*gamma(B*phi,-1.)*Z(Q,llambda)**2
    return (a+b) / (a+b/llambda**2)

# R.H.S of eq. S93
def theta_tilde(f_exc,B,phi,llambda,Q,Z):
    n = f_exc*gamma_prime(B,1.0)*Z(Q,1.)+\
        (1.-f_exc)*phi*gamma_prime(B*phi,-1.0)*Z(Q,llambda)
    a = f_exc*gamma(B,1.)*Z(Q,1.)**2
    b = (1.-f_exc)*gamma(B*phi,-1.)*Z(Q,llambda)**2
    return - n / sqrt(2.*(a+b))

# R.H.S of eq. S91
def C(f_exc,B,phi,llambda,Q,Z):
    return f_exc*gamma(B,1.0)*Z(Q,1.0)**2 + \
           (1.-f_exc)*gamma(B*phi,-1.0)*Z(Q,llambda)**2

@vectorize
def get_OP_kappa_in(kappa_in,p_out,f_exc,phi,llambda):
    '''
    Find order parameters of maximal kappa_in solutions

    ## Paramaters

    kappa_in - input noise robustness (margin)
    p_out - fraction of pattarns with +1 label
    f_exc - fraction of excitatory input afferents
    phi - CV_exc/CV_inh (Eq. S21)
    llambda - sigma_inh/sigma_exc (Eq. S20)

    ## Returns the following order parameters:

    Q, Delta, C, theta_tilde, B, alpha
    '''
    Z=Z_kappa_in(kappa_in,p_out)
    #first we check for the ubalanced solution (see Eqs. S118 and S119):
    Q_sp_equation = lambda Q: Q-Q_eq(f_exc,0.0,phi,llambda,Q,Z)
    #Q = scipy.optimize.brentq(Q_sp_equation,0.1,10.)
    Q = scipy.optimize.newton(Q_sp_equation,1.)
    tt = theta_tilde(f_exc,0.0,phi,llambda,Q,Z)
    #if theta_tilde is greater then zero the solution is unbalanced and B=0
    if tt>0.0:
        B=0.0
    else:
    #Solution is balanced (theta_tilde=0) we need to solve for
    #both Q and B.
        # initial value of Q
        Q0=Q
        # B as a function of Q is the solution to S129
        B_of_Q = lambda QQ: scipy.optimize.brentq(\
                    lambda B: theta_tilde(f_exc,B,phi,llambda,QQ,Z),\
                     -10.,10.)
        # Eq. S128
        Q_sp_equation = lambda QQ: QQ-Q_eq(f_exc,B_of_Q(QQ),phi,llambda,QQ,Z)
        Q = scipy.optimize.newton(Q_sp_equation,Q0)
        # Eq. S129
        B=B_of_Q(Q)
        # for completness we compute theta_tilde (should be zero)
        tt=theta_tilde(f_exc,B,phi,llambda,Q,Z)

    # Aux. var
    k = kappa_in/sqrt(Q)
    # Eq. S127
    D = Delta(k,p_out)
    # Eq. S130
    c = C(f_exc,B,phi,llambda,Q,Z)
    # Eq. S96
    alpha = 2*c*alpha_unconst(D,k,p_out)
    return Q,D,c,tt,B,alpha

@vectorize
def get_OP_kappa_out_balanced(kappa_out,Gamma,p_out,\
                              f_exc,phi,llambda):
    '''
    Find order parameters of maximal kappa_out solutions
    for alpha<alpha_b (balanced solutions)

    ## Paramaters

    kappa_out - output noise robustness
    Gamma - maximal norm of weights
    p_out - fraction of pattarns with +1 label
    f_exc - fraction of excitatory input afferents
    phi - CV_exc/CV_inh (Eq. S21)
    llambda - sigma_inh/sigma_exc (Eq. S20)

    ## Returns the following order parameters:

    Q, Delta, C, theta_tilde, B, alpha

    For these solutions theta_tilde=0
    '''
    Z=Z_kappa_in(kappa_out/Gamma,p_out)
    #Solution is balanced (theta_tilde=0) we need to solve for
    #both Q and B (Equations are the same as the balanced solution of
    #maximal kappa_in solutions).
    B_of_Q = lambda QQ: scipy.optimize.brentq(\
                lambda B: theta_tilde(f_exc,B,phi,llambda,QQ,Z),\
                      -10.,10.)
    Q_sp_equation = lambda QQ: QQ-Q_eq(f_exc,B_of_Q(QQ),phi,llambda,QQ,Z)
    Q = scipy.optimize.newton(Q_sp_equation,1.0)
    B=B_of_Q(Q)
    k = kappa_out/Gamma/sqrt(Q)
    D = Delta(k,p_out)
    c = C(f_exc,B,phi,llambda,Q,Z)
    tt = theta_tilde(f_exc,B,phi,llambda,Q,Z)
    alpha = 2*c*alpha_unconst(D,k,p_out)
    return Q,D,c,tt,B,alpha

#Eq. S137 (TYPO: In this equation f is supposed to be p_out)
def kappa_eq(kappa,B,theta_tilde,c,f_exc,phi,p_out):
    tt=theta_tilde
    k=tt*kappa
    D=Delta(k,p_out)

    return B/sqrt(2*c)+\
           kappa*alpha_unconst(D,k,p_out)*2*\
           (1.-p_out)*(H_prime(D-k)+(D-k)*H(D-k))

# Eq. S138
def C_of_B(B,f_exc,phi):
    a=f_exc*gamma(B,1.0)+(1-f_exc)*gamma(B*phi,-1.0)
    b=f_exc*B*gamma_prime(B,1.0)+(1-f_exc)*B*phi*gamma_prime(B*phi,-1.0)
    return a*(1-b/(2*a))**2

@vectorize
def get_OP_kappa_out_unbalanced(B,p_out,f_exc,phi,llambda):
    '''
    Find order parameters of maximal kappa_out solutions
    for alpha_b<alpha<alpha_c (balanced solutions)

    ## Paramaters

    B - solve for a given OP B. B is zero at capacity and positive below
    p_out - fraction of pattarns with +1 label
    f_exc - fraction of excitatory input afferents
    phi - CV_exc/CV_inh (Eq. S21)
    llambda - sigma_inh/sigma_exc (Eq. S20)

    ## Returns the following order parameters:

    Q, Delta, C, theta_tilde, kappa_0, alpha

    kappa_0 is given in eq. S136
    '''
    # in this case Z does not appear in the SP equations so we define
    # a dummy function that equals 1.
    Z=lambda Q,llambda: 1.0
    # calculat theta_tilde and C for the given B
    tt=theta_tilde(f_exc,B,phi,1.0,1.0,Z)
    c=C_of_B(B,f_exc,phi)
    # find the scaled output robustness
    kappa_0 = scipy.optimize.newton(\
            lambda k: kappa_eq(k,B,tt,c,f_exc,phi,p_out),0.)

    # calculate Delta and alpha
    k=kappa_0*tt
    D = Delta(k,p_out)
    alpha = 2*c*alpha_unconst(D,k,p_out)
    # Finally for completness calculate Q (eq. S132)
    Q = Q_eq(f_exc,B,phi,llambda,1.0,Z)
    return Q,D,c,tt,kappa_0,alpha

def get_alpha_c_and_alpha_b(p_out,f_exc,phi):
    kappa = 0.
    Gamma = 1.
    llambda = 2. # dummy values. alpha_c and alpha_b do not
                       # depend on Gamma or llambda
    Q,D,c,tt,B,alpha_c = get_OP_kappa_in(kappa,p_out,f_exc,phi,llambda)
    Q,D,c,tt,B,alpha_b = get_OP_kappa_out_balanced(kappa,Gamma,p_out,\
                                                   f_exc,phi,llambda)
    return alpha_c,alpha_b


def B_star(f_exc,phi):
    ''' returns the maximal value of B for unbalanced max kappa_out
    soutions.'''
    return scipy.optimize.brentq(\
            lambda B: theta_tilde(f_exc,B,phi,1.0,1.0,lambda Q,llambda: 1.0),\
            0,10.)
