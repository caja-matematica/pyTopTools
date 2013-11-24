import numpy as np
import matplotlib.pyplot as plt
import cPickle as pkl


def normalize(arr, imin=0, imax=1, dmin=None, dmax=None):
    """
    Normalize 'arr', in-place. (Stolen from stack
    overload. Surprised numpy doesn't have a built-in normalize
    function.)

    (imin, imax) -- desired range of normalization

    dmin and dmax -- used if the array does not include all of the
    values. For example, birth time may not include the minimum and
    maximum values. In this case, 0 and max_height are passed to the
    function.
    """
    arr = arr.astype( float )
    if dmin is None:
        dmin = 0 #arr.min()
    if dmax is None:
        dmax = arr.max()
    arr -= dmin
    arr *= (imax - imin) # general range
    arr /= (dmax - dmin)
    arr += imin
    return arr

def find_max( frames ):
    """
    frames : list of persistence diagrams.
    """
    arr = np.vstack( frames )
    return arr.max()

def extract_robust_birth( frame, gen=1 ):
    """
    Diagrams from Perseus should already be ordered by lifespan. We do
    this in case they are not.
    """
    g = gen + 1 # acct for inf generator
    lifespans = np.diff( frame )
    L = lifespans.T[0]
    # argsort give indices of sorted array, from lowest to
    # highest (take only largest 'gen_num' lifespans; final
    # [::-1] reverses order). so the 'infinite' generator
    # lifespan is last. using -g gets us (gen_num+1) slots
    # back. When all is done, the birth time of the infinite
    # generator is at index 0.
    sort_idx = L.argsort()[-g:][::-1]

    # grab the birth coordinate
    births = list( frame[sort_idx,0] )
    return births[1:]


def main( num_gens=1, step=10 ):
    """
    """
    prefix = '/ima/imausr1/jberwald/data/rbc/'
    newfile = 'new_robust.pkl'
    oldfile = 'old_robust.pkl'
    newmax = 'new_cell_maxes.pkl'
    oldmax = 'old_cell_maxes.pkl'

    # load persistence diagrams (cells and frames)
    with open( prefix + newfile ) as fh:
        new = pkl.load( fh )
    with open( prefix + oldfile ) as fh:
        old = pkl.load( fh )

    # load cell maxes
    with open( prefix + newmax ) as fh:
        new_maxes = pkl.load( fh )
    with open( prefix + oldmax ) as fh:
        old_maxes = pkl.load( fh )

    # First normalize birth/death times then extract top N robust
    # generators
    new_var = []
    old_var = []
    for cell, frames in old.items():
        cellmax = None # old_maxes[cell]
        for frame in frames.values():
            nframe = normalize( frame, dmax=cellmax )
            # this returns birth times of robust, non-infinite generators
            robust = extract_robust_birth( nframe, gen=num_gens )
            # finally, compute the variance for this frame
            old_var.append( np.var(robust) )

    for cell, frames in new.items():
        cellmax = None # new_maxes[cell]
        for frame in frames.values():
            nframe = normalize( frame, dmax=cellmax )
            # this returns birth times of robust, non-infinite generators
            robust = extract_robust_birth( nframe, gen=num_gens )
            # finally, compute the variance for this frame
            new_var.append( np.var(robust) )

    # this just creates a prettier graph
    np.random.shuffle( new_var )
    np.random.shuffle( old_var )

    # plot the variance over every [step]'th frame
    fig = plt.figure()
    ax = fig.gca()
    ax.plot( new_var[::step], 'bo', ms=6, label='New variance' )
    ax.plot( old_var[::step], 'r.', ms=6, label='Old variance' )
    xmax = max( len(old_var[::step]), len(new_var[::step]) )
    ax.set_xlim( right=xmax )
    ax.legend()
    
    #labels
    ax.set_xlabel( 'Frame', fontsize=16 )
    ax.set_ylabel( 'Variance of frame', fontsize=16 )

    fig.savefig( './figs_26cells/birththreshold_variance.pdf' )
    fig.show()

    # make box plots of the mean/variance of the variability over all
    # birth times
    new_mean = np.mean( new_var )
    old_mean = np.mean( old_var )
    
    fig2 = plt.figure()
    
    
    return new_var, old_var


if __name__ == "__main__":
    
    new,old = main( 6, step=10 )
