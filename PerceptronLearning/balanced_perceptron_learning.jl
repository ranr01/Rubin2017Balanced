# Julia 0.6.1 code
BLAS.set_num_threads(8)

function perceptronCycle!(w,g,X,y,eps,sigma_in,sigma_out,rho)
    P=length(X)
    err=0
    for ind = randperm(P)
		eta = sqrt(w'*w*sigma_in^2 + sigma_out^2)*randn()
        h = w'*X[ind]-1 + eta
        if y[ind]*h[1]<0.
            err+=1
            w[:] = (1.0-eps) * w + rho*y[ind]*X[ind]
            w[:] = ((g.*w).>0.).*w
        else
			w[:] = (1.0-eps) * w
		end
    end
    return err
end

function sim(N,P,sigma_w_init,sigma_in,sigma_out,
             rho,eps,f_exc,f_in,p_out,T,stat_interval)
    # g is vector of types of neurons
    g=ones(N)
    g[Int(N*f_exc):end]=-1
    # input patterns
    X=[map(Float64,rand(N).>f_in) for μ=1:P]
    # desired output
    y=2*map(Float64,rand(P).>p_out)-1
    # initial weights
    w=sigma_w_init*randexp(N).*g

    TT = Int(T/stat_interval)
    iter=Array{Float64}(TT)
    err=Array{Float64}(TT)
    I=Array{Float64}(TT)
    IB=Array{Float64}(TT)
    w_norm=Array{Float64}(TT)
    kappa_in=Array{Float64}(TT)
    kappa_out=Array{Float64}(TT)

    for t=1:T
        cycle_error = perceptronCycle!(w,g,X,y,eps,
						   			 sigma_out,sigma_in,rho)/P
        if t % stat_interval == 0
            #Gathering stats
            h=Array{Float64}([(w'*x)[1] for x=X])
            hbar=mean(h)
            gw=g.*w
            C=Array{Float64}([(gw'*x)[1] for x=X])
            kappa_mu=y.*(h-1)
            kout=minimum(kappa_mu)
            wn=norm(w)
            ind = Int(t/stat_interval)
            # saving stats
            err[ind] = cycle_error
            kappa_out[ind] = kout
            kappa_in[ind] = kout/wn
            I[ind] = hbar/(sqrt(N)*std(h))
            IB[ind] = hbar/mean(C)
            w_norm[ind] = wn
            iter[ind] = t
        end
    end

    return iter,err,I,IB,w_norm,kappa_in,kappa_out
end

##
P=900
N=3000

sigma_w_init = 1e-3
sigma_out = 0.15
sigma_in = 0.1
rho = 0.1/N
eps = 5e-7

f_in = 0.5
p_out = 0.5

f_exc = 0.8

T = 10000
stat_interval = 10
TT  = Int(T/stat_interval)
n_sim = 2

err1=Array{Float64}(n_sim,TT)
I1=Array{Float64}(n_sim,TT)
w_norm1=Array{Float64}(n_sim,TT)
IB1=Array{Float64}(n_sim,TT)
kappa_in1=Array{Float64}(n_sim,TT)
kappa_out1=Array{Float64}(n_sim,TT)

for n=1:n_sim
    print(n,"\n")
    t,e,i,ib,wn,kin,kout=sim(N,P,sigma_w_init,
                             sigma_in,sigma_out,
							 rho,eps,
                             f_exc,f_in,p_out,
                             T,stat_interval)
    err1[n,:]=e
    I1[n,:]=i
    IB1[n,:]=ib
    w_norm1[n,:]=wn
    kappa_in1[n,:]=kin
    kappa_out1[n,:]=kout
end
##
using Plots
pyplot()
##
p1 = plot(mean(err1,1)',xscale=:log10,label="Error")
p2 = plot(mean(IB1,1)',xscale=:log10,ylim=(0,0.4),label="IB")
p3 = plot(mean(kappa_in1,1)',xscale=:log10,ylim=(0,0.4),label="κ_in")
p4 = plot(mean(kappa_out1,1)',xscale=:log10,ylim=(0,0.4),label="κ_out")
plot(p1,p2,p3,p4,layout=(2,2))
