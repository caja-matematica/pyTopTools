"""
compute_hopf_diagram_dists.py

Author: Jesse Berwald

Opened: August, 2013

This code computes the distances between persistence diagrams for
realizations of a noisy Hopf system. The trial ('run') and window size
('window') are fixed. The routine then compares successive persistence
diagrams at t and t+window.

"""
import numpy as np
import matplotlib.pyplot as plt
import sys, time

# local stuff
from pyTopTools.catastrophes import saddle_node as S
from pyTopTools import timeseries

if len( sys.argv ) < 3:
    print "Meant to be used in batch format, must provide a 'run number'"
    sys.exit(1)

run = int( sys.argv[1] )
window_size = int( sys.argv[2] )
print "run:", run
print "window:", window_size

# for overlapping windows
#step_size = 5

# simple overlappers
step_size = int( window_size / 2 )

# path to realization 'run'
# prefix = '/sciclone/data10/jberwald/climate_attractors/hopfsub/hopfsub_trial'
# # output path for Perseus and distance data
# dist_prefix = '/sciclone/data10/jberwald/climate_attractors/hopfsub_dist_step'+\
#     str(step_size) + '/hopfsub_' 

# for non-overlapping windows
prefix = '/sciclone/data10/jberwald/climate_attractors/hopfsub_data_persStep50/hopfsub_trial'
# output path for Perseus and distance data
dist_prefix = '/sciclone/data10/jberwald/climate_attractors/hopfsub_dist_persStep_50' +\
    '/hopfsub_' 

# for testing
# prefix = './hopfsub_trial'
# dist_prefix = '../data/hopfsub_'

# --> If we are not passing in both run and window size for batch
# processing use the sequence of sizes below. Set of range of window
# sizes appropriate for the 2D system. main issue: huge simplices
# created before the Hopf bifurcation when points are 'tight'. Window
# sizes refer to indices, not 'time' on the t axis window_sizes =
# np.arange( 50, 350, 50 )


#for step_size in window_sizes:
# all distances at this window size
all_distances = []

# read in data from hopfsub_trail{run}.txt 
data = np.loadtxt( prefix + str(run) + '.txt', delimiter=',' )

tvec = data[:,0] # time vector
xy = data[:,1:]  # x,y coords

# for testing
#xy = xy[11000:13000]

tmax = len( xy )

# set up some limits and containers
w0 = int( 0.1 * tmax )
w1 = int( 0.9 * tmax )
previous_window = None
distances = []
nx = []

# for testing
# all_windows = []

# for each realization, move the (non-overlapping) window
# along the time series
while w0 < w1:
    #print "w0:", w0
    left = int( w0 )
    right = int( w0 + window_size )

    # slice the xy vector to get a window
    window = np.array( xy[ left : right ] )

    # compute persistence diagrams. Use WindowND for dim > 1
    current_window = timeseries.WindowND( window, diagram_dim=1 )

    w = left
    nx.append( w )
    persin = dist_prefix + 'w'+ str( w ) + '_step'+\
             str( step_size ) + '_trial'+ str( run ) + '.txt'

    # stepsize and nsteps are fairly snesitive wrt to diagram output
    # as well as computation time! make sure that nsteps is adequate
    # to capture the persistence homological structures.
    current_window.convert2perseus( persin, 
                                    stepsize=0.01, nsteps=50 )
    # make the diagram
    t0 = time.time()
    current_window.compute_persistence( persin )
    print "time to compute PD: ", time.time() - t0

    # if not the first window, compute d(D1,D2)
    if previous_window is not None:
        t0 = time.time()
        dist = current_window.compute_wasserstein_distance( previous_window.perspath )
        print "time to compute Wasserstein distance: ", time.time() - t0
        distances.append( float(dist) )
        
    # update
    previous_window = current_window
    w0 += step_size

all_distances.append( distances )

# make (n realizations) x (windows) array 
distarr = np.asarray( all_distances )

# save an array of  distances for this step_size
np.savetxt( dist_prefix + 'distances_window'+str( window_size ) +'_trial'+str( run )+'.txt', 
            distarr )

# mean_dist = distarr.mean( axis=0 )
# np.savetxt( prefix + 'mean_distances_step'+str( step_size ) + '.txt', mean_dist )



