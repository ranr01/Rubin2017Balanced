import numpy as np
import cvxopt

from types import FunctionType

def interactive(f):
    """decorator for making functions appear as interactively defined.
    This results in the function being linked to the user_ns as globals()
    instead of the module globals().
    """

    # build new FunctionType, so it can have the right globals
    # interactive functions never have closures, that's kind of the point
    if isinstance(f, FunctionType):
        mainmod = __import__('__main__')
        f = FunctionType(f.__code__, mainmod.__dict__,
            f.__name__, f.__defaults__,
        )
    # associate with __main__ for uncanning
    f.__module__ = '__main__'
    return f


def random_patterns_network(N,f_ex,f_inh,g_ex,alpha):
    '''
    Generate P=alpha*N random memory states.
    Each X_i^mu is i.i.d according to:
    X_ex ~ B(f_ex), X_inh ~ B(f_inh).
    The network has N*g_ex excitatory and N*(1-g_ex) inhibitory neurons.

    Returns:
    * X - N x P matrix of 0,1
    * g - X x 1 vector of +-1 according to the neuron's type (E or I)
    '''
    P=int(alpha*N)  # number of memory states
    N_ex=int(g_ex*N) # number of E neurons
    N_inh=N-N_ex # number of I neurons
    X_ex=(np.random.rand(N_ex,P)< f_ex).astype(np.float) # memory patterns of E
    X_inh=(np.random.rand(N_inh,P) < f_inh).astype(np.float) # memory patterns of I
    X=np.vstack((X_ex,X_inh)) # memory pattens
    g=np.ones((N,1)) # type of neurons
    g[N_ex:]=-1.
    return X,g

def learn_network_aff(X,g,i,sim_func):
    '''
    Given desired memory patterns X and neuron types g, find the weigth vector
    of neuron i, w_i, using the fucntions sim_func

    Returns

    * w -  neuron weights assuming threshold = 1
    * converged to solution - boolean. wether sim_func converged to a weight
      implementing all memory states as fixed points.
    '''
    y=2*X[i,:]-1 # desired label of neuron i accross memory states
    X=np.vstack((X[0:i,:],X[i+1:,:])) # memory states removing neuron i
    g=np.vstack((g[0:i,:],g[i+1:,:])) # neuron type's removing neuron i
    #print(X.shape,g.shape,y.shape)
    # Finding weights using sim_func
    cvxopt.solvers.options['show_progress']=False
    w,theta,tau,sol,converged_to_solution=sim_func(X,y,g,Gamma=10.)
    cvxopt.solvers.options['show_progress']=True
    return w/theta, converged_to_solution

@interactive
def learn_network(sim_func,N,lbview):
    '''
    Learns a recurrent network.
    This is meant to be executed on remote engins (Using Ipython Parallel) and
    is desigend to minimize data transfer. Thus it assumes that the memory states,
    X, and neuron types vector are defined globally for each engine
    (this way they only need to be pushed once). In addition each engine must
    also define the sim_func.

    Returns:
    w - a list of N weight vectors (of size N-1) implementing the memory states
    all_conerged - True if sim_func converged to a valid solution for all neurons
    memories_are_fixed_points - True if all memory states are fixed points of the dynamics

    In principle all_converged should be always equal to memories_are_fixed_points.
    If they are not this might indicate a bug.
    '''
    amr=lbview.map(lambda i: learn_network_aff(X,g,i,\
                    sim_func),\
               range(N))
    amr.wait_interactive()
    w,c=zip(*amr)
    all_conerged = np.array(c).all()
    memories_are_fixed_points = array([((dot(w[i].T,\
                    np.vstack((X[0:i,:],X[i+1:,:])))-1.>0)==X[i,:]).all()\
                    for i in range(N)]).all()
    return w, all_conerged, memories_are_fixed_points


class stochastic_network(object):
    '''
    A class for a recurrent network with stochastic dynamics
    '''
    def __init__(self, J):
        '''
            Construct a recurrent network from a weight matrix.
            J is a assumed to be an iterable of length N with weight vector
            of size N-1 each.
        '''
        self.N=len(J[0])+1 # Network size
        # Constructing an N x N connection matrix with J_ii = 0.
        self.J=[np.vstack((j[0:i],[[0]],j[i:])).T \
                for i,j in enumerate(J)]
        self.sigma=1.0 # Gaussian noise magnitude
        self.X=np.zeros((self.N,1)) # initial state of the network.
    def cycle(self):
        '''
            Updates the state of all neurons sequnetialy in random order.
        '''
        for i in np.random.permutation(self.N):
            # calculating the local field. (threshold is 1 for all neurons)
            h = np.dot(self.J[i],self.X)-1
            # updating the state of neuron i
            self.X[i] = float(h+np.random.randn()*self.sigma>0.)

class T0_network(object):
    '''
        A class of a recurrent network with no output noise (zero temperature)
    '''
    def __init__(self, J):
        self.N=len(J[0])+1 # Network size
        # Constructing an N x N connection matrix with J_ii = 0.
        self.J=[np.vstack((j[0:i],[[0]],j[i:])).T \
                for i,j in enumerate(J)]
        self.X=np.zeros((self.N,1)) # initial state of the network.
    def cycle(self):
        '''
            Updates the state of all neurons sequnetialy in random order.
        '''
        for i in np.random.permutation(self.N):
            # calculating the local field. (threshold is 1 for all neurons)
            h = np.dot(self.J[i],self.X)-1
            # updating the state of neuron i
            self.X[i] = float(h>0)

def test_stochastic_dynamics(sigma,J,X_mu,T=500,err_thr=0.15):
    '''
        Tests the stability of a memory state against output noise with SD
        sigma.

        Parameters:

        * sigma - SD of output noise
        * J - Connection matrix. Assumed to be an iterable of length N with
            weight vectors of size N-1 each
        * X_mu - Memory state. N vector of 0,1

        Returns:

        * e - the final error of the network
        * last_t - the time in wich the simulation terminated. For stable patterns
          we expect last_t = T
    '''
    net = stochastic_network(J) #create network
    net.sigma = sigma # set noise level
    N = X_mu.shape[0]
    # initialize network to memory state
    net.X[:] = X_mu
    last_t = T # time of termination of simulation
    for t in range(T):
        # update all neurons
        net.cycle()
        # calculate error
        e = np.sum(np.abs(net.X-X_mu))/N
        if e > err_thr: # no point to continue simulation if to far from memory
            last_t = t+1
            break
    return e,last_t

def noisy_pattern(X,g,f_ex,f_inh,p):
    '''
        Creates a noisy version of memory patterns X with distortion level p
        conserving the mean activity levels of E (f_ex) and I (f_inh) neurons.

        Parameters:
        * X - 0,1 memory patterns (N x P)
        * g - N vector of neuorn types (1 -> E , -1 -> I)
        * f_ex - mean activity level of E neurons
        * f_inh - mean activity level of I neurons
        * p - distortion level of memory patterns
    '''
    N,P=X.shape
    g_mat=g[:,np.zeros(P,dtype=int)]
    P_flip=\
    (g_mat==1)*((X==1)*p+(X==0)*p*f_ex/(1-f_ex))+\
    (g_mat==-1)*((X==1)*p+(X==0)*p*f_inh/(1-f_inh))
    return np.abs(X-(np.random.rand(N,P)<P_flip))

def test_basin(p,J,g,f_ex,f_inh,X_mu,err_thr=0.15):
    '''
        Tests the convergence to a memory state against distortion
        of initial network state.

        Parameters:

        * p - distortion level of memory patterns
        * J - Connection matrix. Assumed to be an iterable of length N with
            weight vectors of size N-1 each
        * g - N vector of neuorn types (1 -> E , -1 -> I)
        * f_ex - mean activity level of E neurons
        * f_inh - mean activity level of I neurons
        * X_mu - Memory state. N vector of 0,1

        Returns:

        * e0 - the initial error of the network
        * e - the final error of the network (after a maximum of 50 network cycles)
    '''
    net = T0_network(J)
    # make sure X_mu is an Nx1 vector
    X_mu = X_mu.reshape((-1,1))
    N = X_mu.shape[0]

    # set initial state to distorted memory pattern
    net.X[:]=noisy_pattern(X_mu,g,f_ex,f_inh,p)
    # calculate initial error
    e0=np.mean(np.abs(net.X-X_mu))
    # run network for at most 50 network cycles.
    for t in range(50):
        net.cycle()
        e=np.mean(np.abs(net.X-X_mu))
        if e==0.0 or e>err_thr: # if network diverged stop
            break
    return e0,e
