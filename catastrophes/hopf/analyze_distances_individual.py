import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
from scipy.stats import linregress

# local
from pyTopTools.catastrophes import saddle_node as S
from pyTopTools import timeseries as T

# run a new stochastic simulation or not
do_sim = False
plot_lr = False
#prefix = '/ima/imausr1/jberwald/data/climate/hopf_step5/hopfsub_distances_window'
prefix = '/ima/imausr1/jberwald/data/climate/hopfsub_dist_persStep50/hopfsub_distances_window'
trial_prefix = '/ima/imausr1/jberwald/data/climate/hopfsub_data_persStep50/hopfsub_trial'
#steps = range( 50, 350, 50 )
#steps = range( 50, 300, 50 )
#steps = [100,150,200]
steps = range( 100, 225, 25 )
trials = range( 1, 31 ) # 30 or 50 trials, depending on experiment

# from Marian's matlab code
dt = 0.05

# window size
window = 150

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

# load all of the trials
hopf_trials = {}
for t in trials[:4]:
    hopf_trials[t] = np.loadtxt( trial_prefix +str(t) + '.txt', delimiter=',' )
    
#
# PLOT SHIT
#
# the gridded subplots
if 0:
    fig, axarr = plt.subplots(2, 1, figsize=(20,8), frameon=False )
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None,
                        wspace=None, hspace=None)
    #plt.axis( 'off' )

    # left plot, full time series 
    ax0 = axarr[0]
    tvec = np.array( hopf_trials[1][:,0] )
    xvec = np.array( hopf_trials[1][:,1] )
    yvec = np.array( hopf_trials[1][:,2] )

    tmax = len( tvec )

    # plot the t,x time series
    # crop the time a bit
    # nx = np.linspace( 0.1*tmax, 0.9*tmax, len(xvec) )
    # nx = np.asarray( nx, dtype=int )
    # tx = tvec[nx]
    ax0.plot( tvec, xvec, lw=2 )
    ax0.set_ylabel( r'$x(t)$', fontsize=16 )
    ax0.set_title( '(Noisy) Hopf bifurcation', fontsize=18 )

    # moving_window = 30
    # for i,s in enumerate( steps[2:] ):
    #     avg = avgs[s]
    #     #nx = np.arange( 0.1*tmax, 0.9*tmax, 5 )[1:-1]
    #     ma = T.moving_average( avg, moving_window )
    #     nx = np.linspace( 0.1*tmax, 0.9*tmax, len(ma) )
    #     #ma = np.log10( ma )
    #     ax0.plot( dt * nx, ma, lw=2, label=r'$w$: '+str(s) )

    ax1 = axarr[1]

    # first distance vec at window=150 lines up with hopf_trila[1]
    dvec = np.array( data[window][0] )
    dstd = np.ones_like( dvec ) * dvec.std() 
    zvec = np.zeros_like( dvec )
    nx = np.linspace( 0.1*tmax, 0.9*tmax, len(dvec) )
    nx = np.asarray( nx, dtype=int )
    tx = tvec[nx]
    #avg = np.log10( avg )
    ax1.plot( tx, dvec, marker='.', lw=2, label=r'window size ($w$): '+str(s) )  
    ax1.fill_between( tx, dstd, facecolor='g', alpha=0.5 )

    if 0:
        divider = make_axes_locatable( ax1 )
        # sub_ax = fig.add_axes([0.6, 0.6, .2, .2], axisbg='y')
        # vals, bins = np.histogram( dvec, 20 )
        # sub_ax.hlines( bins[:-1], [0], vals, linewidth=3 )
        axHisty = divider.append_axes( "right", size=1.2, pad=0.1, sharey=ax1 )
        vals, bins = np.histogram( dvec, 20 )
        axHisty.hlines( bins[:-1], [0], vals, linewidth=10 )
        #axHisty.plot( bins[:-1], vals, '^', ms=10 )
        #axHisty.set_xticklabels( [0, 125, 250], ['0', '125', '250'] )
        axHisty.set_frame_on( False )
        axHisty.set_xticks( [] )
        axHisty.set_yticks( [] )

        # append blank axes
        divider2 = make_axes_locatable( ax0 )
        blank_ax = divider2.append_axes( "right", size=1.2, pad=0.1, sharey=ax0 )
        blank_ax.set_frame_on( False )
        blank_ax.set_xticks( [] )
        blank_ax.set_yticks( [] )

    ax1.set_xlabel( r'window $k$ ($w='+str(window)+')$', fontsize=16 )
    ax1.set_ylabel( r'$W_2$', fontsize=16 )
    ax1.set_xlim( 0, 1000 )
    fig.show()


    # #ax2.set_xlim( left_lim, right_lim +right_buffer )
    # ax2.set_xlim( 0, right_lim +right_buffer )
    # #ax2.set_ylim( ybottom, ytop )
    # ax2.set_xlabel( 'time', fontsize=label_size )
    # ax2.set_ylabel( r'Wasserstein 2-norm $D_i = d_2(P_i^w, P_{i+1}^w)$'+'\n(log scale)', 
    #                 fontsize=label_size )
    # ax2.legend( loc='upper left', fontsize='x-small' )

if 1:
    fig = plt.figure()
    
    # full time series 
    ax0 = fig.add_subplot( 1,1,1 )
    tvec = np.array( hopf_trials[1][:,0] )
    xvec = np.array( hopf_trials[1][:,1] )
    yvec = np.array( hopf_trials[1][:,2] )

    tmax = len( tvec )

    # moving_window = 30
    # for i,s in enumerate( steps[2:] ):
    #     avg = avgs[s]
    #     #nx = np.arange( 0.1*tmax, 0.9*tmax, 5 )[1:-1]
    #     ma = T.moving_average( avg, moving_window )
    #     nx = np.linspace( 0.1*tmax, 0.9*tmax, len(ma) )
    #     #ma = np.log10( ma )
    #     ax0.plot( dt * nx, ma, lw=2, label=r'$w$: '+str(s) 

    # first distance vec at window=150 lines up with hopf_trila[1]
    dvec = np.array( data[window][0] )
    dstd = np.ones_like( dvec ) * dvec.std() 
    zvec = np.zeros_like( dvec )
    nx = np.linspace( 0.1*tmax, 0.9*tmax, len(dvec) )
    nx = np.asarray( nx, dtype=int )
    tx = tvec[nx]
    #avg = np.log10( avg )
    # ax0.plot( tx, dvec, marker='.', lw=2, label=r'window size ($w$): '+str(s) )  
    # ax0.fill_between( tx, dstd, facecolor='g', alpha=0.5 )

    # plot the t,x time series
    # crop the time a bit
    # nx = np.linspace( 0.1*tmax, 0.9*tmax, len(xvec) )
    # nx = np.asarray( nx, dtype=int )
    # tx = tvec[nx]
    scale = 1 #0.5 * dvec.max() / xvec.max()
    line0 = ax0.plot( tvec[::10], scale * xvec[::10], lw=1, label=r'$f(x)$' )
    ax0.set_xlabel( 'time', fontsize=16 )
    ax0.set_ylabel( r'$x(t)$', fontsize=16 )
    ax0.set_title( '(Noisy) Hopf bifurcation', fontsize=18 )
    #ax0.legend( loc='upper left' )

    ax1 = fig.add_subplot( 1,1,1, sharex=ax0, frameon=False )
    line1 = ax1.plot( tx, dvec, 'g', marker='.', lw=3, label=r'$D_{i+1}-D_i$' )  
    ax1.yaxis.tick_right()
    ax1.yaxis.set_label_position("right")
    ax1.set_ylabel( '$D_{i+1}-D_i$' )
    line2 = ax1.fill_between( tx, dstd, facecolor='g', alpha=0.5, label='Std Dev' )

#    ax1.legend( loc='upper left' )
    # plt.figlegend( (line0, line1, line2), 
    #                (r'$f(x)$',r'$D_{i+1}-D_i$', 'Std Dev'), 
    #                'upper left' )

    plt.show()
