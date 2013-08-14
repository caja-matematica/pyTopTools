"""
Compute topological persistence of correlated structures in climate model data

author : Jesse Berwald

opened : July 17, 2013
"""
import numpy as np
import matplotlib.pyplot as plt
from pyTopTools import pyImage

prefix = 'deg'
#prefix = 'bet'
#prefix = 'clus'
#prefix = 'clos'

suffix = '-Temp-H0.2-L0.82-TS60-20sin.dat'
fpath = '/Users/jberwald/Dropbox/Work/Climate/IdeaLab/Data/Temp_H0.2_L0.82_TS60-20sin/'+\
        prefix+'-Temp-H0.2-L0.82-TS60-20sin.dat'


# suffix = '-RBC_ST10_40_TS10'
# fpath = '/Users/jberwald/Dropbox/Work/Climate/IdeaLab/Data/RBC_ST/'+\
#         prefix + suffix + '.dat'


# name for file containing data converted to Perseus format
persin = './data/'+prefix + suffix + '_pers.txt'

#  prefix for Perseus output; appended with '_*.txt', where * is the
#  homology level
persout =  './data/'+prefix + suffix + '_pers'

# homology dimension to analyze
dim = 1

# to get around NaN issue
eps = 1e-15

# scale to blow blow floating points up to integers (see
# write_perseus())
scale = 100

#def persistence( fpath, persin, persout ):
try:
    data = np.load( fpath )
except IOError:
    data = np.loadtxt( fpath )

#data = np.log( data + eps )

# init object to work with Perseus cutting off first row is only for
# the special network data related to heat diffusion
P = pyImage.PerseusImage( data[1:] )

# write data in special perseus format. Assume we have an ndarray
P.write_perseus( persin, dtype='cubtop', scale=scale)

# dtype='cubtop' (cubical toplex, eg., non sparse data);
# dtype='scubtop' (sparse cubical toplex)
P.run_perseus( persout, dtype='cubtop' )

# PLOT SHIT
fig = plt.figure()
ax = fig.gca()
im = ax.imshow( data[1:] / scale, interpolation='nearest' )
fig.colorbar( im )
fig.show()
fig.savefig('/Users/jberwald/Dropbox/Work/Climate/IdeaLab/Data/figures/'+\
            prefix+suffix+'_intensity.png')

# dim refers to the homology dimension we want to analyze
fig0 = P.draw_persdia( 0, scale, ms=0.2 )
ax0 = fig0.gca()
ax0.set_title(r'Connected Components -- $H_0$')
plt.draw()
fig0.savefig('/Users/jberwald/Dropbox/Work/Climate/IdeaLab/Data/figures/'+\
             prefix+ suffix + '_H0.png')

fig1 = P.draw_persdia( 1, scale, ms=0.2 )
ax1 = fig1.gca()
ax1.set_title(r'Bounded regions -- $H_1$')
plt.draw()
fig1.savefig('/Users/jberwald/Dropbox/Work/Climate/IdeaLab/Data/figures/'+
             prefix + suffix + '_H1.png')
