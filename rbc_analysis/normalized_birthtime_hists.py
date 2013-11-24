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

def find_max( frames ):
    """
    frames : list of persistence diagrams.
    """
    arr = np.vstack( frames )
    return arr.max()

def plot_first_robust_hist( new_births, old_births, nbins=100 ):
    """
    """
    fig = plt.figure()
    ax1 = fig.gca()
    
    ax1.hist( new_births, bins=nbins, color='b', edgecolor='none', 
              alpha=0.8, label='First new robust generator', histtype='stepfilled',
              zorder=2 )
    ax1.hist( old_births, bins=nbins, color='r', edgecolor='none', 
              alpha=0.7, label='First old robust generator', histtype='stepfilled',
              hatch='.', zorder=2 )
    
    new_mean = np.mean( new_births )
    old_mean = np.mean( old_births )
    ax1.vlines( new_mean, 0, ax1.get_ylim()[1], linestyle='dashed', 
                label='New mean', zorder=1 )
    ax1.vlines( old_mean, 0, ax1.get_ylim()[1], linestyle='dashed', 
                label='Old mean', zorder=1 )
    
    # labels
    ax1.set_xlabel( 'Normalized birth threshold', fontsize=16 )
    ax1.set_ylabel( 'Threshold frequency', fontsize=16 )
    
    ax1.legend()
    plt.show()

    fig.savefig( './figs_26cells/birththresholds_top_overlay.pdf' )


def main( num_gens=1, nbins=100 ):
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

    # load cell avgs
    with open( prefix + newmax ) as fh:
        new_maxes = pkl.load( fh )
    with open( prefix + oldmax ) as fh:
        old_maxes = pkl.load( fh )


    # First N normalized birth/death times
    new_normed = []
    old_normed = []
    for cell, frames in old.items():
        avg = old_maxes[cell]
        for frame in frames.values():
            nframe = normalize( frame, dmax=avg )
            # this returns birth times of robust, non-infinite generators
            robust = extract_robust_birth( nframe, gen=num_gens )
            old_normed.append( robust )

    for cell, frames in new.items():
        avg = new_maxes[cell]
        for frame in frames.values():
            nframe = normalize( frame, dmax=avg )
            # this returns birth times of robust, non-infinite generators
            robust = extract_robust_birth( nframe, gen=num_gens )
            new_normed.append( robust )

    # combine all birth thresholds to get a cumulative list -->
    # histogram
    new_births = []
    old_births = []
    for x in new_normed:
        new_births += [ xi for xi in x ]
    for y in old_normed:
        old_births += [yi for yi in y ]    

    # extract 1st, 2nd, 3rd generators  and plot histograms
    fig1 = plt.figure()
    ax1 = fig1.gca()

    # now plot complete histograms
    ax1.hist( new_births, bins=nbins, color='b', edgecolor='none', 
              alpha=0.7, label='New, top '+str(num_gens)+'  generators', 
              histtype='stepfilled', zorder=2 )
    ax1.hist( old_births, bins=nbins, color='r', edgecolor='none', 
              alpha=0.7, label='Old, top '+str(num_gens)+' generators', 
              histtype='stepfilled', hatch='.', zorder=2 )

    # plot new and old birth threshold means *behind* histograms
    new_mean = np.mean( new_births )
    old_mean = np.mean( old_births )
    ax1.vlines( new_mean, 0, ax1.get_ylim()[1], linestyle='dashed', 
                label='New mean', zorder=1 )
    ax1.vlines( old_mean, 0, ax1.get_ylim()[1], linestyle='dashdot', 
                label='Old mean', zorder=1 )

    ax1.set_xlabel( 'Normalized birth threshold', fontsize=16 )
    ax1.set_ylabel( 'Threshold frequency', fontsize=16 )
    ax1.legend()    
    plt.show()

    fig1.savefig( './figs_26cells/all_birththresholds_top'+str(num_gens)+'.pdf' )
        
    return new_normed, old_normed

if __name__ == "__main__":
    
    new,old = main( 3, nbins=100 )

    new_births = []
    old_births = []
    # extract the first non-infinite generator
    for x,y in zip( new, old ):
        try:
            new_births.append( x[0] )
        except IndexError:
            pass
        try:
            old_births.append( y[0] )
        except IndexError:
            pass
    plot_first_robust_hist( new_births, old_births )
