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
from pyTopTools import timeseries


def euler_maruyama( f, g, dt, x, lam, sigma ):
    """
    dX = f( x, lam )dt + g( x, sigma )dW
    
    args = [ x (spatial), lam (parameter), noise (sigma) ]
    """
    return x + dt * f( x, lam ) + sqrt( dt ) * gauss( 0, 1 ) * g( x, sigma )

def integrate( f, g, h, x0, tmax, dt, lam0=0, sigma=0 ):
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

def sde_integrator( f, g, h, tmax=10000, xinit=-2, dt=0.01, 
                    sigma=0.2, lam0=-4 ):
    """
    #=======================
    # Equation 2.2 (double-well potential with noise)
    #=======================
    f = lambda x, lam : lam + x - 0.005*x**3
    g = lambda x, sigma : sigma # * x (spatial variability)
    h = lambda x, lam : 0.001  #* x 
    
    """    
    #transient = int( 0.5 * tmax/dt )
    tvec, xvec, lamvec = integrate( f, g, h, xinit, tmax, dt,
                                    lam0=lam0, sigma=sigma )
    return tvec, xvec, lamvec


def plot_timeseries( nt, nx, color='b', lw=1 ):
    """
    """
    fig = plt.figure()
    ax = fig.gca()
    ax.plot( nt, nx, 'b-', c=color, lw=lw, ms=2 )
    plt.show()
    return fig


if __name__ == "__main__":

    f = lambda x, lam : lam + x - 0.005*x**3
    g = lambda x, sigma : sigma # * x (spatial variability)
    h = lambda x, lam : 0.001  #* x 

    tmax = 10000
    xinit = -2
    dt = 0.01
    sigma = 0.2
    lam0 = -4
    transient = int( 0.5 * tmax/dt )

    tvec, xvec, lamvec = integrate( f, g, h, xinit, tmax, dt,
                                    lam0=lam0, sigma=sigma )

    fig = plot_timeseries( tvec[transient:], xvec[transient:] )
    # fig.savefig( 'general_bif_dt'+str(dt)+'.png' )

    # fig2 = plt.figure()
    # ax2 = fig2.gca()
    # drift = []
    # vf = np.vectorize( f )
    # nx = np.linspace( -15, 15, 1000 )
    # for l in lamvec[::10000]:
    #     drift.append( vf( nx, l ) )
    # for d in drift:
    #     ax2.plot( d )
    # fig2.show()

    # far from bifurcation
    i0 = int( 2020 / dt )
    i1 = int( 2050 / dt )
    a = np.array( xvec[ i0 : i1 ] )
    a -= a.mean()

    # right next to bifurcation
    j0 = int( 9370 / dt )
    j1 = int( 9400 / dt )
    b = np.array( xvec[ j0 : j1 ] )
    b -= b.mean()

    # fig3 = plt.figure()
    # ax3 = fig3.gca()
    # ax3.hist( a, 20 )
    # ax3.set_xlabel( r'state ($x(t)$)', fontsize=14 )
    # ax3.set_ylabel( r'Number of occurences', fontsize=14 )
    # fig3.savefig('figures/state_distribution_far_from_bif.png', transparent=True )
    # fig3.show()

    # fig4 = plt.figure()
    # ax4 = fig4.gca()
    # ax4.hist( b, 20 )
    # ax4.set_xlabel( r'state ($x(t)$)', fontsize=14 )
    # ax4.set_ylabel( r'Number of occurences', fontsize=14 )
    # fig4.savefig('figures/state_distribution_near_bif.png', transparent=True )
    # fig4.show()

    fig5 = plt.figure()
    ax5 = fig5.gca()
    ax5.plot( a, 'b-', lw=2, label='Far from transition' )
    ax5.plot( b, 'r--', lw=2, label='Near transition' )
    ax5.set_xlabel( r'$t$ (time step in window)', fontsize=14 )
    ax5.set_ylabel( r'$x(t)$', fontsize=14 )
    ax5.legend()
    fig5.savefig('figures/state_near_far.png', transparent=True )
    fig5.show()

    # # compute persistence diagrams
    # window1 = timeseries.Window( a )
    # window2 = timeseries.Window( b )

    # vrfile = './data/sadde_'
    # window1.convert2perseus( vrfile + 'far.txt', stepsize=0.0005, nsteps=50 )
    # window2.convert2perseus( vrfile + 'near.txt', stepsize=0.0005, nsteps=50 ) 

    # # make the diagrams
    # window1.compute_persistence( vrfile + 'far.txt' )
    # window2.compute_persistence( vrfile + 'near.txt' )

    # # plot the diagrams
    # f1 = window1.draw_diagram( vrfile + 'far_0.txt', scale=True, marker_scale=0.05 )
    # f2 = window2.draw_diagram( vrfile + 'near_0.txt', scale=True, marker_scale=0.05 )
    # f2.savefig('figures/saddle_persdia_far.png', transparent=True)
    # f2.savefig('figures/saddle_persdia_near.png', transparent=True)
