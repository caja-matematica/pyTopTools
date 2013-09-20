import numpy as np
import matplotlib.pyplot as plt
import sys

# local stuff
from pyTopTools.catastrophes import saddle_node as S
from pyTopTools import timeseries

dataname = './deltaO18_2200_2400.txt'
vrdata = './vr/deltaO18_'

middle = 200

# "near" the bifurcation
data1 = np.loadtxt( './deltaO18_2350_2450.txt', delimiter=',' )  
# "far" from the bifurcation
data2 = np.loadtxt( './deltaO18_2450_2550.txt', delimiter=',' )
# "further"
data3 = np.loadtxt( './deltaO18_2550_2650.txt', delimiter=',' )

#
# first window
#
window1 = np.array( data1 )
window1 = window1[::-1]
window1 -= window1.mean() # detrend
W1 = timeseries.Window( window1 )
persin1 = vrdata + 'w1.txt'
W1.convert2perseus( persin1,
                    stepsize=0.0005, nsteps=50 )
# make the diagram
W1.compute_persistence( persin1 )

#
# second window 
#
window2 = np.array( data2 )
window2 = window2[::-1]
window2 -= window2.mean() # detrend
W2 = timeseries.Window( window2 )
persin2 = vrdata + 'w2.txt' 
W2.convert2perseus( persin2,
                    stepsize=0.0005, nsteps=50 )
# make the diagram
W2.compute_persistence( persin2 )

#
# third window 
#
window3 = np.array( data3 )
window3 = window3[::-1]
window3 -= window3.mean() # detrend
W3 = timeseries.Window( window3 )
persin3 = vrdata + 'w3.txt' 
W3.convert2perseus( persin3,
                    stepsize=0.0005, nsteps=50 )
# make the diagram
W3.compute_persistence( persin3 )


# analyze stuff
dist12 = W1.compute_wasserstein_distance( W2.perspath )
dist23 = W2.compute_wasserstein_distance( W3.perspath )

dists = [ int(dist12), int(dist23) ]

print "distance 1 to 2: ", dist12
print "distance 2 to 3: ", dist23

figs = []
for i,w in enumerate( [W1,W2,W3] ):
    fig = w.draw_diagram() 
    fig.savefig( '../figures/deltaO18_w'+str(i)+'.pdf', transparent=True )
    figs.append( fig )
    
F = plt.figure()
for subf in figs:
    F.axes.append( subf.gca() )
    
