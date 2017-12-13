import cvxopt
import numpy as np

def sign_constrained_perceptron_max_kappa_out(X,y,g,Gamma=1.,Lambda=1e5,\
                                                  external_input=None):
    '''Finds the maximal $\kappa_\mathrm{out}$ solution for a sign constrained
        Percptron, with |w|<=Gamma.

        ## Paramaters:

        X - Input patterns (N x P)
        y - Patterns' labels (P x 1)
        g - Sign of weights (N x 1, of +1 for excitatory and -1 for inhibitory)
        Gamma - maximal norm of solution's weight vector
        Lambda - Regularization parameter for feasibitilty variable (Should be >>1)

        ## Returns

        w - Perceptron weights
        theta - Perceptron threshold
        tau - Regularization variable (tau=0 if solution found and tau>0 if no
        solution exists)
        sol - Full output dictionary from the CVXOPT solver
        converged_to_solution - True if solution found (based on actual classification)

        NOTE: w and theta are the so called canonical weights. The weigths in units
        of threshold are given by threshold * w / theta.
        In terms of w and theta $\kappa_\mathrm{out}$ is given by 1/theta and
        $\kappa_\mathrm{in}$ is given by 1/|w|.
    '''
    N,P =X.shape
    y = np.array(y).reshape((P,1))
    g = np.array(g).reshape((N,1))

    if external_input is None:
        Theta_mu = np.ones((P,1))
    else:
        Theta_mu = -np.array(external_input).reshape((P,1))+1.

    A = np.hstack([X.T*y[:,np.zeros(N,int)]*g[:,np.zeros(P,int)].T, \
                 -y*Theta_mu, np.ones((P,1))])
    beta = np.ones((P,1))
    a = np.zeros((N+2,1))
    a[N] = 1.
    a[N+1] = Lambda
    #We need to solve min(a^Tx) subject to:
    # Ax>=beta
    # x>=0
    # |x|<Gamma*theta

    #cvxopt solves:
    # min(c^Tx) subjet to
    # Gx+s=h
    # s>=0

    #second order cone:
    # s0=Gamma*x[N]
    # i=1...N si=x[i-1]
    # s0>=||s||
    # So matrix is N+1XN+2
    G_0 = np.zeros((N+1,N+2))
    G_0[0,-2] = -Gamma
    for i in range(1,N+1):
        G_0[i,i-1] = -1

    c = cvxopt.matrix(a)
    G = cvxopt.matrix(np.vstack(\
                [np.diag(np.vstack([-np.ones((N,1)),[[-1.],[-1]]]).flatten()),\
                 -A,\
                 G_0])\
                )

    h = cvxopt.matrix(np.vstack([np.zeros((N+2,1)),-beta,np.zeros((N+1,1))]))
    dims = {'l':N+2+P, 'q': [N+1], 's':[]}

    # Solving the linear program
    sol = cvxopt.solvers.conelp(c,G,h,dims)

    # extracting solution
    w = g*np.array(sol['x'][:N])
    theta = sol['x'][N]
    tau = sol['x'][N+1]

    # testing the solution
    converged_to_solution = (y.T*(np.dot(w.T,X)-theta*Theta_mu.T)>=1.).all()

    if not converged_to_solution:
        print("Did not find solution. tau={}".format(tau))
    return w,theta,tau,sol,converged_to_solution

def sign_constrained_perceptron_max_kappa_in(X,y,g,Gamma=1.0,Lambda=1e5,\
                                             external_input=None):
    '''Finds the maximal $\kappa_\mathrm{in}$ solution for a sign constrained
        Percptron, with |w|<=Gamma.

        ## Paramaters:

        X - Input patterns (N x P)
        y - Patterns' labels (P x 1)
        g - Sign of weights (N x 1, of +1 for excitatory and -1 for inhibitory)
        Gamma - maximal norm of solution's weight vector
        Lambda - Regularization parameter for feasibitilty variable (Should be >>1)

        ## Returns

        w - Perceptron weights
        theta - Perceptron threshold
        tau - Regularization variable (tau=0 if solution found and tau>0 if no
        solution exists)
        sol - Full output dictionary from the CVXOPT solver
        converged_to_solution - True if solution found (based on actual classification)

        NOTE: w and theta are the so called canonical weights. The weigths in units
        of threshold are given by threshold * w / theta.
        In terms of w and theta $\kappa_\mathrm{out}$ is given by 1/theta and
        $\kappa_\mathrm{in}$ is given by 1/|w|.
    '''
    N,P = X.shape
    y = np.array(y).reshape((P,1))
    g = np.array(g).reshape((N,1))

    if external_input is None:
        Theta_mu = np.ones((P,1))
    else:
        Theta_mu = -np.array(external_input).reshape((P,1))+1.

    A = np.hstack([X.T*y[:,np.zeros(N,int)]*g[:,np.zeros(P,int)].T, \
                -y*Theta_mu, np.ones((P,1))])
    beta = np.ones((P,1))
    a = np.zeros((N+2,1))
    a[N+1] = Lambda
    Q = np.eye(N+2)
    Q[N,N] = 0.0
    Q[N+1,N+1] = 0.0

    # We need to solve min(1/2x^TQx+a^Tx) subject to:
    # Ax>=beta
    # x>=0
    # |x|<Gamma*theta

    #second order cone:
    # s0=Gamma*x[N]
    # i=1...N si=x[i-1]
    # s0>=||s||
    # So matrix is N+1XN+2
    G_0 = np.zeros((N+1,N+2))
    G_0[0,-2] = -Gamma
    for i in range(1,N+1):
        G_0[i,i-1] = -1

    #cvxopt qp solves:
    # min(1/2x^TPx+q^Tx) subjet to
    # Gx+s=h
    # S>=0
    q = cvxopt.matrix(a)
    P_opt = cvxopt.matrix(Q)
    G = cvxopt.matrix(np.vstack(\
                [np.diag(-np.ones(N+2)),\
                 -A,\
                 G_0]))
    h = cvxopt.matrix(np.vstack([np.zeros((N+2,1)),-beta,np.zeros((N+1,1))]))
    dims = {'l':N+2+P, 'q': [N+1], 's':[]}
    #print(dims)

    # solving the quadratic program
    sol = cvxopt.solvers.coneqp(P_opt,q,G,h,dims)

    # extracting solution
    w = g*np.array(sol['x'][:N])
    theta = sol['x'][N]
    tau = sol['x'][N+1]

    # testing the solution
    converged_to_solution = (y.T*(np.dot(w.T,X)-theta*Theta_mu.T)>=1.).all()
    if not converged_to_solution:
        print("Did not find solution. tau={}".format(tau))
    return w,theta,tau,sol,converged_to_solution
