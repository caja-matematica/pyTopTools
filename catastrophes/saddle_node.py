'''
saddle_node.py

Author: Jesse Berwald

Opened: April 11, 2013

Created data for a 1D saddle node bifurcation, with noise in phase
space and drift in the parameter.
'''
import numpy as np
from random import gauss
from math import sqrt
import matplotlib.pyplot as plt


def euler_maruyama( f, g, dt, x, lam, sigma ):
    """
    dX = f( x, lam )dt + g( x, sigma )dW
    
    args = [ x (spatial), lam (parameter), noise (sigma) ]
    """
    return x + dt * f( x, lam ) + sqrt( dt ) * gauss( 0, 1 ) * g( x, sigma )

def integrate( f, g, h, x0, nsteps, dt, lam0=0, sigma=0 ):
    """
    Integrate the stochastic differential equation (Langevin eqn) 

    dX = f( x, lam )dt + g( x, sigma )dW
    d\lambda = h( x, lam )dt

    where f is a deterministic function, g is the stochastic portion,
    x is a variable in phase space, and lambda is a parameter.
    """
    nx = [ x0 ]
    nt = [ 0. ]
    lams = [ lam0 ]

    # loop optimizations
    sqdt = sqrt( dt )
    trange = np.linspace( 0, tmax, tmax/dt ) 
    for i in trange:
        # integrate according to Euler-Maruyama
        step = euler_maruyama( f, g, dt, nx[-1], lams[-1], sigma )
        
        # evolve lambda deterministically to get drift
        lstep = lams[-1] + dt * h( nx[-1], lams[-1] )
        
        nx.append( step )
        lams.append( lstep )
        nt.append( nt[-1] + dt )

    return nt, nx, lams

def plot_timeseries( nt, nx, color='b', lw=1 ):
    """
    """
    fig = plt.figure()
    ax = fig.gca()
    ax.plot( nt, nx, '-o', c=color, lw=lw, ms=2 )
    plt.show()
    return fig


if __name__ == "__main__":

    f = lambda x, lam : lam + x - x**3
    g = lambda x, sigma : sigma # * x (spatial variability)
    h = lambda x, lam : 0.05  #* x 

    tmax = 150
    xinit = -2
    dt = 0.01
    sigma = 0.1
    lam0 = -6

    tvec, xvec, lamvec = integrate( f, g, h, xinit, tmax, dt,
                                    lam0=lam0, sigma=sigma )



    fig = plot_timeseries( tvec, xvec )
    #fig.savefig( 'general_bif_dt'+str(dt)+'.png' )
