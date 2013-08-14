import numpy as np
import matplotlib.pyplot as plt
import cPickle as pkl
from collections import defaultdict
import time

# local
from pyTopTools import chomp_tools as C


def check_frame( frame, cone=None, bound=1, plot_dia=False, 
                 cell=None, frame_num=None):
    """
    """
    if cone is None:
        cone = np.median(frame[:,0]) 

    lifespans = np.diff( frame ).ravel()
    w = np.where( lifespans > bound )
    longlived = frame[w]
    #longlived = frame
        
    # DIAGRAM WITH CONE
    if plot_dia:
        plot_diagram_cone( frame, cone, bound, cell=cell, num=frame_num )

    # born before
    w0 = np.where( longlived[:,0] < cone )[0]
    # dies after
    w1 = np.where( longlived[:,1] > cone )[0]
    
    robust_idx = np.intersect1d( w0, w1 )
    robust = longlived[ robust_idx ]
    
    return robust

def plot_diagram_cone( dia, cutoff, bound, cell=None, num=None,
                       save_prefix='/sciclone/home04/jberwald/github/local/caja-matematica/pyTopTools/rbc_analysis/figs_26cells/' ):
    """
    """
    line = lambda x : x + 100
    nx = np.linspace( 0, dia.max() )
    ny = line( nx )
    plt.figure()
    diag = [ 0, dia.max() ]
    plt.plot( diag, diag, 'g-', lw=2 )
    plt.plot( nx, ny, 'r-', lw=1 )
    plt.plot( dia[:,0], dia[:,1], 'bo' )
    plt.hlines( cutoff, 0, cutoff, linestyle='dashed' )
    plt.vlines( cutoff, cutoff, dia.max(), linestyle='dashed' )
    plt.show()
    # if cell is not None:
    #     plt.savefig( save_prefix + cell + '_betti_cone_f'+str(num)+'.pdf' )
    
def get_robust_betti( data, cone='median', bound=1, max_scale=0.3,
                      new_maxes=None, old_maxes=None,
                      plot_dia=False ):
    """Track threshold betti number corresponding to robust features

    data : filename, eg. new_robust.pkl or old_robust.pkl

    *_maxes : file containing the maximum pixel value over frames or
    entire cells

    singel_cell : Should we only compute for a single cell
    """
    try:
        print " Reading file..."
        # data is a file name
        with open( data ) as fh:
            cells = pkl.load( fh )
        cell_vals = cells.values()
    except TypeError:
        # in this case we assume that data contain a list of a single
        # cell's generators in each frame
        cell_vals = [ data ]

    # read in maxes if we want to normalize values
    if old_maxes is not None:
        with open( norm_old ) as fh:
            maxes = pkl.load( fh )
    elif new_maxes is not None:
        with open( norm_new ) as fh:
            maxes = pkl.load( fh )
    else:
        maxes = None

    # hold all the feature information (birth time, death time, or
    # longest lifespan)
    # if longest_lifespan:
    features = defaultdict( list )

    # for each cell, stack the generators into one array
    for cell, cell_frames in cells.items():
        print "cell", cell

        # list of robust threshold arrays in this frame
        cell_vals = cell_frames.values()

        # compute the median birth value over all robust generators
        if cone is 'median':
            all_vals = np.vstack( cell_vals )
            med = np.median( all_vals[:,0] )

        # used to cutoff relative to max value
        if maxes is not None:
            cell_max = maxes[ cell ]

        # each frame contains (n x 2) array of (birth,death) thresholds
        for i, frame in enumerate( cell_vals ):
            robust = check_frame( frame, cone=med, bound=bound, 
                                  plot_dia=plot_dia, cell=cell, frame_num=i )
            features[cell].append( len(robust) )

    return features

def compute_difference_intervals( betti_ts ):
    """Compute intervals between changes in robust betti numbers.
    """
    d = np.diff( betti_ts )

    # print 'D mean', d.mean()
    # dw = np.asarray( d, dtype=np.float )
    # dw -= dw.mean()

    w = np.where( d != 0 )[0]
    dw = np.diff( w )

    print "dw mean", dw.mean()
    dw = np.asarray( dw, dtype=np.float )
    dw -= dw.mean()
    return dw

def plot_robust_interbetti_mse( fname, eps ):
    """
    """
    with open( fname ) as fh:
        B = pkl.load( fh )
    
    eps = str(eps)

    new_mse = []
    old_mse = []

    for cell, cell_betti in B.items():
        print cell
        if cell.startswith( 'new' ):
            color = 'b'
        else:
            color = 'r'

        dw = compute_difference_intervals( cell_betti )
        
        np.savetxt( 'data/'+cell+'_betti_interbeat_eps'+str(eps)+'.txt', dw )
        C.run_mse( 'data/'+cell+'_betti_interbeat_eps'+str(eps)+'.txt' )

        # give plenty of time to write result to disk
        time.sleep( 1 )

        mse = C.mse_converter( 'data/'+cell+'_betti_interbeat_eps'+str(eps)+'.mse' )
        
        if cell.startswith('new'):
            new_mse.append( mse[:,1] )
        else:
            old_mse.append( mse[:,1] )

        # ax = plt.subplot2grid( (2,1), (0,0) )
        # # interbeat plot
        # ax.plot( dw, color+'-' )
        # ax.set_xlabel( 'Change' )
        # ax.set_ylabel( r'Time interval between robust $\beta_1$ changes' )
        # ax.set_title( 'MSE of intervals between changes in robust generators -- '+cell )
        # ax = plt.subplot2grid( (2,1), (1,0) )
        # # MSE plot
        # ax.plot( mse[:,0], mse[:,1], 'm-', lw=2 )
        # ax.set_ylabel( 'Sample entropy') 
        # ax.set_xlabel( 'Resolution' )
        # ymin, ymax = ax.get_ylim()
        # ax.set_ylim( 0, ymax )
        # plt.savefig( 'figs_26cells/'+cell+'_betti_interbeat_MSE_'+eps+'.pdf' )

    # grab the last x-axis values
    nx = mse[:,0]
    newarr = np.vstack( new_mse )
    oldarr = np.vstack( old_mse )
    new_mean = newarr.mean( axis=0 )
    old_mean = oldarr.mean( axis=0 )
    
    fig = plt.figure()
    ax = fig.gca()
    ax.plot( nx, new_mean, 'b-', lw=2, label="Mean MSE, young" )
    ax.plot( nx, old_mean, 'r-', lw=2, label="Mean MSE, old" )   
    ax.set_title( 'Mean MSE of intervals between changes in robust generators' )
    ax.set_ylabel( 'Sample Entropy' )
    ax.set_xlabel( 'Resolution' )
    ax.legend()
    plt.show()
    plt.savefig( './figs_26cells/mean_MSE_interbeat_betti'+eps+'.pdf' )

def plot_robust_betti_mse( fname, eps ):
    """
    fname : path to dictionary { cellname : # robust betti per frame }
    """
    with open( fname ) as fh:
        B = pkl.load( fh )
    
    new_mse = []
    old_mse = []
    eps = str(eps)
    for cell, cell_betti in B.items():
        print cell
        if cell.startswith( 'new' ):
            color = 'b'
        else:
            color = 'r'
        np.savetxt( 'data/'+cell+'_betti_ts.txt', cell_betti )
        C.run_mse( 'data/'+cell+'_betti_ts.txt' )
        time.sleep( 1 )
        mse = C.mse_converter( 'data/'+cell+'_betti_ts.mse' )

        if cell.startswith('new'):
            new_mse.append( mse[:,1] )
        else:
            old_mse.append( mse[:,1] )

        ax = plt.subplot2grid( (2,1), (0,0) )
        ax.plot( cell_betti, color+'-' )
        ax.set_xlabel( 'time' )
        ax.set_ylabel( r'# Robust $\beta_1$' )
        ax.set_title( 'MSE for # of robust generators -- '+cell )
        ax = plt.subplot2grid( (2,1), (1,0) )
        ax.plot( mse[:,0], mse[:,1], 'm-', lw=2 )
        ax.set_ylabel( 'Sample entropy') 
        ax.set_xlabel( 'Resolution' )
        ymin, ymax = ax.get_ylim()
        ax.set_ylim( 0, ymax )
        plt.savefig( './figs_26cells/'+cell+'_robust_betti_MSE_'+eps+'.pdf' )
        
    # grab the last x-axis values
    nx = mse[:,0]
    newarr = np.vstack( new_mse )
    oldarr = np.vstack( old_mse )
    new_mean = newarr.mean( axis=0 )
    old_mean = oldarr.mean( axis=0 )
    
    fig = plt.figure()
    ax = fig.gca()
    ax.plot( nx, new_mean, 'b-', lw=2, label="Mean MSE, young" )
    ax.plot( nx, old_mean, 'r-', lw=2, label="Mean MSE, old" )   
    ax.set_title( 'Mean MSE for number of robust generators' )
    ax.set_ylabel( 'Sample Entropy' )
    ax.set_xlabel( 'Resolution' )
    ax.legend()
    plt.show()
    plt.savefig( './figs_26cells/mean_MSE_robust_betti_'+eps+'.pdf' )

    # new_arr = np.concatenate( new )
    # old_arr = np.concatenate( old )

    # # fignew = plt.figure()
    # # figold = plt.figure()
   
    # plt.figure( figsize=(16,10) )
    # for i,x in enumerate( new[:4] ):
    #     ax = plt.subplot2grid( (4,1), (i,0) )
       
    #     ax.plot( x, 'b-' )
    # plt.figure( figsize=(16,10) )
    # for i,x in enumerate( old[:4] ):
    #     ax = plt.subplot2grid( (4,1), (i,0) )
    #     ax.plot( x, 'r-' )

    # plt.show()


           
if __name__ == "__main__":

    prefix='/sciclone/data10/jberwald/RBC/cells/persout/'

    bnd = 1
    eps = 60
    cellpath = prefix + 'robust_eps'+str(eps)+'.pkl'
    bettis = get_robust_betti( cellpath, bound=bnd, plot_dia=False )
    with open( prefix + 'robust_betti_eps'+str(eps)+'.pkl', 'w' ) as fh:
        pkl.dump( bettis, fh )
    
        
    
        
            
