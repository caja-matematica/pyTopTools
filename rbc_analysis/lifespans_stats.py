import numpy as np
import cPickle as pkl
import matplotlib.pyplot as plt
import time 

from pyTopTools import chomp_tools as C

fdir = '/Users/jberwald/github/local/caja-matematica/pytoptools/rbc_analysis/data/longest_lifespans/'
newfile = 'new_robust_longest_lifespan.pkl'
oldfile = 'old_robust_longest_lifespan.pkl'

with open( fdir + newfile ) as fh:
    new = pkl.load( fh )

with open( fdir + oldfile ) as fh:
    old = pkl.load( fh )

min_new = 5000
for k,c in new.items():
    if len(c) < min_new:
        min_new = len(c)

min_old = 5000
for k,c in old.items():
    if len(c) < min_old:
        min_old = len(c)

min_length = min( min_new, min_old )

new_arr = np.empty( (13, min_length), dtype=int )
old_arr = np.empty( (13, min_length), dtype=int )

for i,v in enumerate( new.values() ):
    new_arr[i] = v[ :min_length ]
for i,v in enumerate( old.values() ):
    old_arr[i] = v[ :min_length ]

new_mse = []
for c in new_arr:
    np.savetxt( 'data/MSE_TMP.txt', c )
    C.run_mse( 'data/MSE_TMP.txt' )
    
    # give plenty of time to write result to disk
    time.sleep( 0.2 )
    mse = C.mse_converter( 'data/MSE_TMP.mse' )
    
    new_mse.append( mse )
    

old_mse = []
for c in old_arr:
    np.savetxt( 'data/MSE_TMP.txt', c )
    C.run_mse( 'data/MSE_TMP.txt' )
    
    # give plenty of time to write result to disk
    time.sleep( 0.2 )
    mse = C.mse_converter( 'data/MSE_TMP.mse' )
    
    old_mse.append( mse )

new_mse_stack = np.vstack( [new_mse[i][:,1] for i in range(len(new_mse))] )
old_mse_stack = np.vstack( [old_mse[i][:,1] for i in range(len(old_mse))] )
new_mean_mse = new_mse_stack.mean( axis=0 )
old_mean_mse = old_mse_stack.mean( axis=0 )

nx = new_mse[0][:,0]

plt.figure()
plt.plot( nx, new_mean_mse, new_mean_mse[1], 'b-' )
plt.plot( nx, old_mean_mse, old_mean_mse[1], 'r-' )



# new_mean = new_arr.mean( axis=0 )
# new_median = np.median( new_arr, axis=0 )

# old_mean = old_arr.mean( axis=0 )
# old_median = np.median( old_arr, axis=0 )

# plt.figure()
# plt.title( 'Robust lifespan means' )
# plt.plot( new_mean, 'b-', lw=2 )
# plt.plot( old_mean, 'r-', lw=2 )

# plt.figure()
# plt.title( 'Robust lifespan median' )
# plt.plot( new_median, 'b-', lw=2 )
# plt.plot( old_median, 'r-', lw=2 )
# plt.show()
