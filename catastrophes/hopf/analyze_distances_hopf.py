import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import linregress

# local
from pyTopTools.catastrophes import saddle_node as S
from pyTopTools import timeseries as T

# run a new stochastic simulation or not
do_sim = False
plot_lr = False
prefix = '/ima/imausr1/jberwald/data/climate/hopf/hopfsub_distances_step'
steps = range( 50, 350, 50 )
trials = range( 1, 51 )

# limits for linear regression and zoom window.
left_lim = 500
#middle = ??

# right hand limit for lin. reg.
right_lim = 700

# additional 
right_buffer = 0 # 1000 - right_lim

# plot settings
label_size = 12 # fontsize

# load all of the data
data = {}
for s in steps:
    st = []
    for t in trials:
        try:
            st.append( np.loadtxt( prefix + str(s) + '_trial'+str(t) + '.txt' ) )
        except:
            continue
    data[s] = st 


# averages for each different window size
avgs = {}
for s,d in data.items():
    v = np.vstack( d )
    avgs[s] = v.mean( axis=0 )

# a static picture of what's going on (bif. around 9400 in this
# realization)
print "Loading data from disk..."
allvec = np.loadtxt( '/ima/imausr1/jberwald/data/climate/hopfsub_example.txt',
                     delimiter=',' )
tvec = allvec[:,0]
xvec = allvec[:,1]
yvec = allvec[:,2]

# set max x-axis value (time)
#tmax = tvec.max()
tmax = len( tvec )

# # linear regression slope to each distance array, for t \in
# # (8000,9000). returns slope, intercept, r_value, p_value, std_err
# if plot_lr:
#     line = lambda slope, intercept, xvals : slope * xvals + intercept
#     linreg1 = {}
#     linreg2 = {}
#     for s in steps:
#         nx = np.arange( 0.1*tmax, 0.9*tmax, s )[:-1]
#         t0 = np.where( nx >= left_lim )[0]
#         t1 = np.where( nx <= middle )[0]
#         t2 = np.where( nx >= middle )[0]
#         t3 = np.where( nx <= right_lim )[0]
#         # two windows for the LR line, with the rh one showing a steep increase.
#         t_int1 = np.intersect1d( t0, t1 )
#         t_int2 = np.intersect1d( t2, t3 )
#         avg = avgs[s]
#         # linear regression line to left and right hand portion
#         m1,b1,r1,p1,std_err1 = linregress( nx[ t_int1 ], avg[ t_int1 ] )
#         m2,b2,r2,p2,std_err2 = linregress( nx[ t_int2 ], avg[ t_int2 ] )
#         print "step:", s
#         print "slope:", m1, m2
#         print ""
#         # linear regression line to  hand portion
#         linreg1[s] = np.array( zip( nx[ t_int1], line( m1, b1, nx[ t_int1 ] ) ) )
#         linreg2[s] = np.array( zip( nx[ t_int2], line( m2, b2, nx[ t_int2 ] ) ) )
    

# plot the little system
#fig = plt.figure(1)
# zoom to a region defined by t \in (left_lim, right_lim)
tvec = np.asarray( tvec )
xvec = np.asarray( xvec )

# left interval for time
tleft = np.where( tvec >= left_lim )[0]
tright = np.where( tvec <= right_lim + right_buffer )[0]
subt = np.intersect1d( tleft, tright )

# # both time intervals for zooming
# subt = np.union1d( subt_1, subt_2 )

# figure with 2x2 subplot structure
fig, axarr = plt.subplots(2, 2, figsize=(20,8) )
# left plot, full time series 
axarr[0,0].plot( tvec, xvec, 'b-' )
#axarr[0,0].set_xlabel( 'time', fontsize=16 )
axarr[0,0].set_ylabel( r'$x(t)$', fontsize=16 )
axarr[0,0].set_title( 'Hopf bifurcation', fontsize=18 )
#
# right plot, uses subt = subt_1 \cup subt_2
# ZOOMED region
#
axarr[0,1].plot( tvec[ subt ], xvec[ subt ], 'b.-', ms=2 )
axarr[0,1].set_xlabel( 'time', fontsize=label_size )

# CHANGE THIS TO ALTER THE BIFURCATION REGION ZOOM IMAGE
# axarr[0,1].set_ylim( -12, 0 )
axarr[0,1].set_title( 'Bifurcation region', fontsize=18 )

####
# plot persistence diagram distance statistics
####
#fig2, axarr2 = plt.subplots(1, 2, figsize=(20,8) ) #plt.figure(2, figsize=(10,8) )
ax2 = axarr[1,0]

# using steps array we add the persdis distance plots in increasing
# step size
for i,s in enumerate( steps[2:] ):
    avg = avgs[s]
    nx = np.linspace( 0.1*tmax, 0.9*tmax, len(avg) )
    nx = np.asarray( nx, dtype=int )
    tx = tvec[nx]
    #avg = np.log10( avg )
    ax2.plot( tx, avg, marker='.', lw=2, label=r'window size ($w$): '+str(s) )  

#ax2.set_xlim( left_lim, right_lim +right_buffer )
ax2.set_xlim( 0, right_lim +right_buffer )
#ax2.set_ylim( ybottom, ytop )
ax2.set_xlabel( 'time', fontsize=label_size )
ax2.set_ylabel( r'Wasserstein 2-norm $D_i = d_2(P_i^w, P_{i+1}^w)$'+'\n(log scale)', 
                fontsize=label_size )
ax2.legend( loc='upper left', fontsize='x-small' )

#
# plot 'derivatives'
#
ax3 = axarr[1,1]
for i,s in enumerate( steps[2:] ):
    delta = np.diff( avgs[s] )
    nx = np.linspace( 0.1*tmax, 0.9*tmax, len(delta) )
    nx = np.asarray( nx, dtype=int )
    tx = tvec[nx]
    ax3.plot( tx, delta, marker='.', lw=2, label=r'$w$: '+str(s) )

ax3.yaxis.tick_right()
ax3.yaxis.set_label_position( 'right' )
ax3.set_xlim( left_lim, right_lim + right_buffer )
#ax3.set_ylim( ybottom, ytop )
ax3.set_xlabel( 'time', fontsize=label_size )
ax3.set_ylabel( r'Derivative: $\Delta_i = D_i - D_{i-1}$',
                fontsize=label_size )
ax3.legend( loc='lower left',fontsize='x-small' )

plt.show()
