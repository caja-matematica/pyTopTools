from pyTopTools.rbc_analysis import birth_times as bt
import matplotlib.pyplot as plt
import cPickle as pkl
import numpy as np


prefix = '/data/jberwald/rbc/timeseries/'

oldname = 'oldgen_avgshift_lb'
newname = 'newgen_avgshift_lb'

new_single = 'newgen_avgshift_lb40.pkl'
old_single = 'oldgen_avgshift_lb40.pkl'

new_cells = [ 'new_50125', 'new_40125', 'new_130125' ]
old_cells = [ 'old_100125', 'old_50125', 'old_90125' ]

# c_new = 'new_50125'
# c_old = 'old_100125'

# options for bt.characteristic_birth_times()
lifespans = True
top = True

gens = [1,2,,3#,4]
val = 40  # lower bound shouldn't matter, so just grab the lb=40 file

ts_new = {}
ts_old = {}

for g in gens:

    print "generator", g
    
    #for val in lb:
    # with open( prefix+ newname +str(val)+'_gen'+str(g)+'_V2.pkl' ) as fh:
    #     new = pkl.load(fh)
    #     ts_new[g] = new

    # with open( prefix+ oldname +str(val)+'_gen'+str(g)+'_V2.pkl' ) as fh:
    #     old = pkl.load(fh)
    #     ts_old[g] = old

    ####################
    ## FOR OVERLAYING SINGLE CELLS
    with open( prefix+ new_single ) as fh:
        # load all birth times for new cells with lb=40
        a = pkl.load(fh)
        for c_new in new_cells:
            sing_new = a[ c_new ]
            bt_new = bt.characteristic_birth_time( sing_new, gen_num=g )
            bt_new = np.array( bt_new )
            bt_new_stand = bt_new = bt_new - bt_new.mean()
            np.savetxt( '/data/jberwald/rbc/birth_time_ts/gen'+str(g)+'_'+c_new+'.txt', bt_new )
            np.savetxt( '/data/jberwald/rbc/birth_time_ts/gen'+str(g)+'_'+c_new+'_mean.txt',
                        bt_new_stand )

    with open( prefix+ old_single ) as fh:
        b = pkl.load(fh)
        for c_old in old_cells:
            sing_old = b[ c_old ]
            bt_old = bt.characteristic_birth_time( sing_old, gen_num=g )
            bt_old = np.array( bt_old )
            bt_old_stand = bt_old - bt_old.mean()
            np.savetxt( '/data/jberwald/rbc/old_gen'+str(g)+'_'+c_old+'.txt', bt_old )
            np.savetxt( '/data/jberwald/rbc/birth_time_ts/gen'+str(g)+'_'+c_old+'_mean.txt',
                        bt_old_stand )

   
    #
    ####################

    new = np.array( new, dtype=int )
    old = np.array( old, dtype=int )

    fig = rh.plot_hist( new, nbins=100 )
    fig = rh.plot_hist( old, nbins=100, fig=fig, color='r' )
    
    fig = rh.plot_hist( bt_new, nbins=100, fig=fig, color='g' )
    fig = rh.plot_hist( bt_old, nbins=100, fig=fig, color='y' )

    # fig = plt.figure()
    ax = fig.gca()
    # ax.set_title( r'Birth times for robust generator #'+str(g)+\
    #                 ' (lifespan threshold = '+str(val)+')', fontsize=16 )
    # ax.set_xlabel( r'Birth threshold (' +c_new+' and '+ c_old +')', fontsize=16 )
    # ax.set_ylabel( r'Number of generators ($g='+str(g)+'$)', fontsize=16 )
    plt.show()
    # fig.savefig( './data/birth_times_gen'+str(g)+'_' + c_old +'_'+c_new + '.png' )
