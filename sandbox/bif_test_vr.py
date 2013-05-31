"""
This module explores bifurcation detection around a saddle node
bifurcation using VR complexes for 'random' points (the y values of a
time series) scattered in 1D.
"""

from pyTopTools import timeseries
from pyTopTools import bottleneck_distance as BD
import matplotlib.pyplot as plt
from numpy import array 

data = timeseries.Timeseries( 'genbif0.txt', 'vr' )

block = 50
stepsize = 0.0001
nsteps = 500

pname = '/Users/jberwald/github/local/caja-matematica/pyTopTools/sandbox/vr_data/genbif_pers'
last_window = None
bdists = []
wdists = []
for t in range( 4000, 8600, block ):
    window = timeseries.Window( data.data, t, t+block )
    window.convert2perseus( pname + str(t)+'.txt', stepsize=stepsize, nsteps=nsteps )
    window.compute_persistence( pname + str(t)+'.txt', pname + str(t), )
    # fig = window.draw_diagram( pname + str(t) + '_0.txt' )
    # ax = fig.gca()
    # ax.set_title( 'genbif -- t0=' +str( t ) +'; block='+str(block) )
    # timeseries.plt.draw()
    # fig.savefig( pname + str(t) + '_0.png' )

    if last_window is not None:
        #d = window.compute_bottleneck_distance( last_window.persdia, engine='py' )
        bdist = window.compute_bottleneck_distance( last_window,
                                                    pname+str(t)+'_0.txt',
                                                    engine='c' )
        wdist = window.compute_wasserstein_distance( pname+str(t)+'_0.txt',
                                                     last_window )
                                                
        wdists.append( (t, wdist) )
        bdists.append( (t, bdist) )
    last_window = pname + str(t) + '_0.txt' 

print "bottleneck distances: ", bdists
print "wasserstein distances: ", wdists

fig2 = plt.figure()
ax = fig2.gca()
b = array( bdists )
b = b.T
w = array( wdists )
w = w.T
ax.plot( b[0], b[1], 'ro' )
ax.plot( w[0], w[1], 'bo-', ms=4 )
fig2.show()
