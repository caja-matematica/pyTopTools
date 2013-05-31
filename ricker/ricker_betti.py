#!/usr/bin/python
"""
Explore the evolution of Betti numbers over time for various fitness
and dispersal values.

Call ricker_beeti.py -h for help menu.
"""
import ricker_patch_model as ricker
import numpy as np
import argparse
import matplotlib.pyplot as plt
import matplotlib.animation as animation

parser = argparse.ArgumentParser(description='Run Ricker coupled patch model on various fitness and dispersal values.')
parser.add_argument( '-r', '--fitness', dest='fitness', type=int,
                   help='A fitness value (integer > 0).', default=5 )
parser.add_argument( '-d', '--dispersal', dest='dispersal', type=float,
                   help='Dispersal value (float in [0,1]).', default=0 )
parser.add_argument( '-m', '--moran', dest='moran', action='store_true',
                     help="Compute Moran's I for spatial correlation. [Default=False]" )
parser.add_argument( '-t', '--tfinal', dest='tfinal', type=int,
                     help="Length of simulation (number of iterations).", default=10 )
parser.add_argument( '-s', '--spinup', dest='spinup', type=int,
                     help='Spin up time (number of transient iterations to discard). Must be < tfinal.',
                     default=0 )
parser.add_argument( '-n', '--patches', dest='patches', type=int,
                     help="Number of patches on one side of matrix. n^2 = size of Ricker matrix.",
                     default=10 )

args = parser.parse_args()

fitness = args.fitness
dispersal = args.dispersal
moran = args.moran
tfinal = args.tfinal
spinup = args.spinup
size = args.patches
prefix = '/data/ricker/betti/'

#ic = 10 * np.random.random( (size,size) )
ic = np.load('n100_ic.npy')

# remove '.' from dispersal so chomp can read the file
if not dispersal==0:
    dd = str( dispersal ).split('.')[1]
else:
    dd = str( dispersal )
betti_name = prefix + 'ricker_r'+str( fitness ) + '_d'+str( dd ) + \
    '_T'+str( tfinal ) + '_n'+str( size )

    #betti_name = None

out = ricker.run( ic, tfinal, fitness, dispersal, spin_up=spinup, save_pop=True,
                  ndim=size, betti_name=betti_name, compute_moran=moran )

b = out['betti']
b0 = b[:,0]
b0 = b0[1::2]
b1 = b[:,1]
b1 = b1[1::2]

nx = np.arange( len(b0) )

fig = plt.figure()
ax = fig.gca()
ax.plot( nx, b0, 'bs-', label=r'$\beta_0$' )
ax.plot( nx, b1, 'ro-', label=r'$\beta_1$' )
ax.set_xlabel('time')
ax.set_ylabel(r'$\beta_0,\quad \beta_1$', fontsize=16)
ax.set_ylim(0,800)
ax.legend( loc='center right' )
fig.show()

# r = out['ricker']
# r.draw_matrix()
# plt.show()


# animate binarized population matrices
# B = out['binary']
# fig2 = plt.figure()
# ims = []
# for x in B:
#     ims.append( (plt.imshow( x, interpolation='nearest' ),) )

# image_animation = animation.ArtistAnimation( fig2, ims, interval=40, repeat_delay=3000 )
# image_animation.save( 'ring_test.mp4', fps=5, codec='ffmpeg' )

#plt.show()
