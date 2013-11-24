import numpy as np
import matplotlib.pyplot as plt

cellpath = '/Users/jberwald/Dropbox/Papers/BioMath/rbc/figures/new11_frame2000.npy'
row = 100

# load the cell frame
a = np.load( cellpath )
a.resize( (203,198) )

# slice
s = a[ row ]

fig = plt.figure()
ax = fig.gca()
ax.plot( s, lw=3, label='Membrane intensity' )
ax.legend()
#ax.set_xticks( [], [] )

# minor dips range from [487,912] and [515,1169] for
# new11_frame2000.npy
ax.hlines( 487, 0, 55, linestyles='dashdot' )
ax.hlines( 912, 0, 35, linestyles='dashdot' )


plt.show()
