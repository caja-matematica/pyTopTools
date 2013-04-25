"""
Module for plotting histograms of generator data.
"""
import numpy
import matplotlib.pyplot as plt
from matplotlib.mlab import stineman_interp


slash = '/'


def plot_hist( ts, nbins=1000, color='b', xlim=None, alpha=0.6, fig=None, **kwargs ):
    """
    ts -- single time series of values to bin.

    nbins -- number of bins to use

    color -- color of the bins

    xlim -- tuple, (xmin,xmax)

    Returns Figure object.
    """
    if not fig:
        fig = plt.figure( figsize=(12,8) )
        fig.patch.set_alpha( 0.0 )
    ax = fig.gca()
    if hasattr( ts, "__array__" ):
        n, bins, patches = ax.hist( ts, bins=nbins, color=color, log=True,
                                    edgecolor='none', alpha=alpha, **kwargs )
        plot_interp = False
    else:
        all_n = []
        for d in ts:
            n, bins, patches = ax.hist( d, bins=nbins, color=color, log=True,
                                        edgecolor='none', alpha=alpha, **kwargs )
            all_n.append( n )

        arr = numpy.asarray( all_n )
        avg = arr.mean( axis=0 )
        err = arr.mean( axis=0 )
        xi, y_interp = interpolate( avg, bins )
        plot_interp = True

    ax.set_xlabel( 'Lifespan', fontsize=20 )
    ax.set_ylabel( 'Number of generators', fontsize=20 )
    ax.tick_params( axis='both', which='major', labelsize=16 )
    #ax.tick_params(axis='both', which='minor', labelsize=8)
    ax.set_ylim( bottom=0.5 ) # accounts for log scale

    if plot_interp:
        ax.plot( xi, y_interp, 'r-', lw=3 )

    if xlim:
        ax.set_xlim( xlim[0], xlim[1] )

    return fig


def interpolate( ts, bins ):
    yp = None
    xi = numpy.linspace( 0, bins[-1], 200 )
    yi = stineman_interp( xi, bins[:-1], ts, yp )
    # interpolate upper and lower error bars to get envelope
    # COMMENTED OUT BELOW
    # upper_yi = stineman_interp( xi, bins[:-1], upper, yp )
    # lower_yi = stineman_interp( xi, bins[:-1], lower, yp )
    return xi, yi
                        
if __name__ == "__main__":

    import cPickle as pkl

    fdir = '/Users/jberwald/github/local/caja-matematica/pyRBC/data/'
    new = fdir + 'new_hist_ts.pkl'
    old = fdir + 'old_hist_ts.pkl'

    with open(new) as fh:
        A = pkl.load(fh)
      
    with open(old) as fh:
        B = pkl.load(fh)

    A.extend(B)
    fig= rh.plot_hist( A, nbins=200 )
    fig= rh.plot_hist( Y, nbins=200 )

