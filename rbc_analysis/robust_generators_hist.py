import birth_times as bt
import rbc_histogram as rh
import matplotlib.pyplot as plt
import cPickle as pkl
import numpy as np


"""
This module is used to create histograms of robust birth times.
"""

def birth_time_hist( prefix, save_prefix):
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
    

def birth_times_hist2( old_path, new_path, save_prefix, gen_num,
                       old_single=None, new_single=None, 
                       old_robustgens_path=None, new_robustgens_path=None, 
                       old_avg_path=None, new_avg_path=None, 
                       death=False, lspans=False, shift=False, 
                       normed=False, # for the concatenated histograms
                       norm_new=None, norm_old=None ):
    """{old,new}_path : path to {old,new}_robust_gen*

    {old,new}_single : name of single old/new cell to overlay
    
    {old,new}_single_path : path to old/new robust generators. gen_num'th
    generator will be extracted using
    birth_times.characteristic_birth_times() function.

    save_prefix : full path to dump images to. will be appended to
    with gen_num.

    {old,new}_avg : path to file containing dictionary of pixel height
    average for each cell

    """
    g = gen_num
    plot_overlay = False

    # container for figure objects
    all_figs = []

    ##########
    # load full path to robust cell dictionaries
    ##########
    with open( new_path ) as fh:
        new = pkl.load(fh)
        
    with open( old_path ) as fh:
        old = pkl.load(fh)

    return new, old

    # if not ( old_single is None or new_single is None ):
    #     ##########
    #     # load single cells
    #     ##########
    #     with open( new_single_path ) as fh:
    #         a = pkl.load(fh)
    #     sing_new = a[ new_single ]
        
    #     with open( old_single_path ) as fh:
    #         b = pkl.load(fh)
    #     sing_old = b[ old_single ]

    ####################
    ## FOR OVERLAYING SINGLE CELLS
    
    # returns the g'th generator, g \in gens. Automatically
    # excludes 'infinite' generator since g \ge 1, and the
    # inf. gen. will sit in position 0 after argsort() in
    # characteristic_birth_times()
    # SINGLE NEW CELL
    #bt_new = bt.characteristic_birth_time( sing_new, g )
    bt_new = bt.characteristic_generator_data_shift( new_robustgens_path, g, single_cell=new_single, 
                                                     birth_times=True, norm_new=norm_new )
    bt_new = np.array( bt_new )

    # subtract the cell avg
    if shift:
        with open( new_avg_path ) as fh:
            new_avg_dict = pkl.load( fh )
        new_avg = new_avg_dict[ new_single ]
        bt_new = np.array( bt_new ) - new_avg

    # SINGLE OLD CELL
    #bt_old = bt.characteristic_birth_time( sing_old, g )
    bt_old = bt.characteristic_generator_data_shift( old_robustgens_path, g, single_cell=old_single, 
                                                     birth_times=True, norm_old=norm_old )
    bt_old = np.array( bt_old )
    # load the avg
    if shift:
        with open( old_avg_path ) as fh:
            old_avg_dict = pkl.load( fh )
        old_avg = old_avg_dict[ old_single ]
        # subtract the cell average
        bt_old = np.array( bt_old ) - old_avg

    plot_overlay = True

    #
    ####################

    if not normed:
        new = np.array( new, dtype=int )
        old = np.array( old, dtype=int )
    else:
        new = np.array( new )
        old = np.array( old )

    new_min = new.min()
    old_min = old.min()
    hist_min = min( new_min, old_min )
    # # shift the bins >= 0
    # new -= hist_min
    # old -= hist_min
    # bt_new -= hist_min
    # bt_old -= hist_min


    # account for single *really large* negative number in old
    # generators
    HUGE = -1e5
    w = np.where( old > HUGE )
    old = old[w]

    # for plotting vertical lines below
    new_mean = new.mean()
    old_mean = old.mean()

    if normed:
        log = False

    # plot the histogram of all birth times for generator g
    fig = rh.plot_hist( old, nbins=100, color='r', hatch='/', log=log, label='Old' )    
    fig = rh.plot_hist( new, nbins=100, fig=fig, log=log, label='Young' )
    
    # plot histogram for single cells (old and new)
    if plot_overlay:
        fig = rh.plot_hist( bt_new, nbins=20, fig=fig, color='g', alpha=0.9, 
                            label='Single (young)', log=log )
        fig = rh.plot_hist( bt_old, nbins=20, fig=fig, color='y', alpha=0.9, 
                            label='Single (old)', log=log )

    # make correct labels
    if death:
        navg_label = 'New death threshold mean'
        oavg_label = 'Old death threshold mean'
        nxlabel = 'Death threshold'
    elif lspans:
        navg_label = 'New lifespans mean'
        oavg_label = 'Old lifespans mean' 
        nxlabel = 'Lifespans'
    else:
        navg_label = 'New birth threshold mean'
        oavg_label = 'Old birth threshold mean'
        nxlabel = 'Birth threshold'
        if normed:
            nxlabel += ' (normalized)'

    ax = fig.gca()
    ymin, ymax = ax.get_ylim()
    # plot some means
    ax.vlines( new_mean, ymin, ymax, linestyle='dashed', label=navg_label )
    ax.vlines( old_mean, ymin, ymax, linestyle='dashdot', label=oavg_label )

    # ax.set_title( r'Birth times for robust generator #'+str(g)+\
    #                 ' (lifespan threshold = '+str(val)+')', fontsize=16 )
     
    #ax.set_xlabel( nxlabel+' (' +new_single+' and '+ old_single +')', fontsize=16 )
    ax.set_xlabel( nxlabel, fontsize=16 )
    if not lspans:
        ax.set_ylabel( r'Number of generators ($g='+str(g)+'$)', fontsize=16 )
    else:
        ax.set_ylabel( r'Number of generators', fontsize=16 )
    ax.legend( loc='upper right' )

    if death:
        plot_type = 'death_times'
    elif lspans:
        plot_type = 'lifespans'
    else:
        plot_type = 'birth_times'
    if normed:
        plot_type += '_normed'
    if shift:
        plot_type += '_shift'
    
    # assume both overlay celsl provided
    if plot_overlay:
        fig.savefig( save_prefix + plot_type + '_gen'+str(g)+'_' + old_single +\
                     '_'+new_single + '.pdf' )
    else:
        fig.savefig( save_prefix + plot_type + '_gen'+str(g)+'.pdf' )
    all_figs.append( fig )
    return all_figs
    
