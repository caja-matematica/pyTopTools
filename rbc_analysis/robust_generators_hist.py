import birth_times as bt
import rbc_histogram as rh
import matplotlib.pyplot as plt
import cPickle as pkl
import numpy as np


"""
This module is used to create histograms of robust birth times.
"""

def birth_time_hist( prefix, save_prefix ):
    #prefix = '/data/jberwald/rbc/timeseries/'

    # append generator number to this filename below...
    oldname = 'oldgen_avgshift_lb'
    newname = 'newgen_avgshift_lb'

    new_single = 'newgen_avgshift_lb40.pkl'
    old_single = 'oldgen_avgshift_lb40.pkl'
    
    # TODO: add more cells to the mix!
    new_cells = [ 'new_50125', 'new_40125', 'new_130125', 'new_140125' ]
    #new_cells = [ 'new_50125' ]
    old_cells = [ 'old_100125', 'old_50125', 'old_90125', 'old_120125' ]
    #old_cells = [ 'old_120125' ]

    gens = [1,2,3] #,4]
    val = 40  # lower bound doesn't matter, so just grab the lb=40 file

    all_figs = []
    for c_new, c_old in zip( new_cells, old_cells ):
        for g in gens:
            # load full histograms
            with open( prefix+ newname +str(val)+'_gen'+str(g)+'_V2.pkl' ) as fh:
                new = pkl.load(fh)

            with open( prefix+ oldname +str(val)+'_gen'+str(g)+'_V2.pkl' ) as fh:
                old = pkl.load(fh)

            ####################
            ## FOR OVERLAYING SINGLE CELLS
            with open( prefix+ new_single ) as fh:
                a = pkl.load(fh)
                sing_new = a[ c_new ]

            with open( prefix+ old_single ) as fh:
                b = pkl.load(fh)
                sing_old = b[ c_old ]

            # returns the g'th generator, g \in gens. Automatically
            # excludes 'infinite' generator since g \ge 1, and the
            # inf. gen. will sit in position 0 after argsort() in
            # characteristic_birth_times()
            bt_new = bt.characteristic_birth_time( sing_new, gen_num=g )
            bt_new = np.array( bt_new )
            bt_old = bt.characteristic_birth_time( sing_old, gen_num=g )
            bt_old = np.array( bt_old )
            #
            ####################

            new = np.array( new, dtype=int )
            old = np.array( old, dtype=int )

            fig = rh.plot_hist( new, nbins=100 )
            fig = rh.plot_hist( old, nbins=100, fig=fig, color='r' )

            fig = rh.plot_hist( bt_new, nbins=100, fig=fig, color='g' )
            fig = rh.plot_hist( bt_old, nbins=100, fig=fig, color='y' )

            ax = fig.gca()
            # ax.set_title( r'Birth times for robust generator #'+str(g)+\
            #                 ' (lifespan threshold = '+str(val)+')', fontsize=16 )
            ax.set_xlabel( r'Birth threshold (' +c_new+' and '+ c_old +')', fontsize=16 )
            ax.set_ylabel( r'Number of generators ($g='+str(g)+'$)', fontsize=16 )
            fig.savefig( save_prefix + 'birth_times_gen'+str(g)+'_' + c_old +'_'+c_new + '.pdf' )
            all_figs.append( fig )
    return all_figs
    
