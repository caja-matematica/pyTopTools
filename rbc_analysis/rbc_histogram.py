"""
Module for plotting histograms of generator data.
"""
import numpy
import matplotlib.pyplot as plt
from matplotlib.mlab import stineman_interp
from rbc_utils import *

slash = '/'


def plot_hist( ts, nbins=1000, color='b', xlim=None, alpha=0.6, fig=None, log=True, **kwargs ):
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
        n, bins, patches = ax.hist( ts, bins=nbins, color=color, log=log,
                                    edgecolor='none', alpha=alpha, **kwargs )
    else:
        all_n = []
        for d in ts:
            n, bins, patches = ax.hist( d, bins=nbins, color=color, log=log,
                                        edgecolor='none', alpha=alpha, **kwargs )
            all_n.append( n )

    ax.set_xlabel( 'Lifespan', fontsize=20 )
    ax.set_ylabel( 'Number of generators', fontsize=20 )
    ax.tick_params( axis='both', which='major', labelsize=16 )
    #ax.tick_params(axis='both', which='minor', labelsize=8)
    ax.set_ylim( bottom=0.5 ) # accounts for log scale

    if xlim is not None:
        ax.set_xlim( xlim[0], xlim[1] )
        
    return fig


def plot_hist_all( ts, nbins=50, transparent=True, norm_it=False, **kwargs ):
    """
    ts -- dictionary with values as lifespan timeseries.

    kwargs -- see pylab hist() function

    Returns interpolation functions as well as histogram triple and figure instance.
    """
    from pylab import log1p
    from matplotlib.mlab import stineman_interp

    #data = ts.values()
    data = ts
    fig = plt.figure( figsize=(12,8) )
    if transparent:
        fig.patch.set_alpha( 0.0 )
    # ax = fig.add_subplot( 121 )
    # ax2 = fig.add_subplot( 122 )
    ax2 = fig.add_subplot( 111 )
    # now plot a single cell's histogram on the first axis
    #ax.hist( data, **kwargs )

    # convert each list of lifespans to 1d arrays
    data = [ numpy.asarray( d ) for d in data ]
    xmax = max( [d.max() for d in data] )
    thebins = numpy.linspace(0, xmax, nbins )

    # add/change the kwargs 'bins' key
    kwargs['bins'] = thebins

    # holds bins counts (n) for each histogram
    all_ny = []
    for d in data:
        n, bins, patches = ax2.hist( d, **kwargs )
        all_ny.append( n )
        # all_bins.append( bins )
    # convert bins counts to a single array to find min
    arr = numpy.array( all_ny ).ravel()
    wy = arr[ numpy.where( arr != 0 )[0] ]
    min_y = min( wy )

    # for plotting average -- these are already log values if
    # log==True in kwargs
    yhist = numpy.array( all_ny, dtype=numpy.float64 )#.ravel()
    avg = yhist.mean( axis=0 )
    err = yhist.std( axis=0 )
    upper = avg + err
    lower = avg - err

    # print "yhist", yhist
    # print "avg", avg
    # print ""

    # print err
    # print ""
    # print lower
    # print upper
    # print ""

    # average value for each histogram bin
    for i, x in enumerate( avg ):
        if x == 0.0:
            avg[i] = 1.0

    # label the axes
    ax2.set_xlabel( 'Generator lifespan', fontsize=20 )
    ax2.set_ylabel( 'Number of generators (log)', fontsize=20 )
    xticks = ax2.get_xticks()
    yticks = ax2.get_yticks()
    ax2.set_xticklabels( [str(int(x)) for x in xticks], fontsize=20 )
    ax2.set_yticklabels( [str(int(y)) for y in yticks], fontsize=20 )

    # now plot the interpolated average distribution here so it is on
    # top of the other stuff
    yp = None
    xi = numpy.linspace( 0, bins[-1],200)
    yi = stineman_interp( xi, bins[:-1], avg, yp )
    # interpolate upper and lower error bars to get envelope
    # COMMENTED OUT BELOW
    upper_yi = stineman_interp( xi, bins[:-1], upper, yp )
    lower_yi = stineman_interp( xi, bins[:-1], lower, yp )
    
    # make sure lower does not go negative since this makes no sense.
    for i,v in enumerate( lower_yi ):
        if v < 1.0:
            lower_yi[i] = 1.0

    # make sure that the plot doesn't get messed up by small values
    # (esp. New cells)
    masked_yi = numpy.ma.masked_less( yi, 1 )

    # plot the interpolation of the avg and the envelope
    ax2.plot( xi, masked_yi, 'r-', lw=3 )
    # ax2.fill_between( xi, lower_yi, upper_yi, #where=masked_yi
    #                   color='r', alpha=0.5, zorder=10 )
    #fig.show()

    if norm_it:
        y_max = masked_yi.max()
        print "y_max", y_max
        masked_yi /= y_max

    return xi, masked_yi, lower_yi, upper_yi, fig, (n, bins, patches)


def plot_all_cell_hists( data, norm_it=True, nbins=100, num_pts=1000, cell_type='new' ):
    """Same as plot_hist_overlay() below, except that along with the mean
    over all cells, the outline of the histogram over each individual
    cell is plotted. Also, the single cell histogram is neglected.

    num_pts : number of points for the interpolation function
    """
    fig = plt.figure()
    ax = fig.gca()

    # lifespans for each cell
    for lifespans in data:
        ny, nx = numpy.histogram( data, bins=nbins, density=True )
        # nx = numpy.arange( 0, num_pts )  
        # fx = numpy.interp( nx, bins[:-1], ny )
  
        # plot the data for this cell
        ax.plot( nx[:-1], ny, lw=2 )

    #plt.ylim( 0.00001, 0.1 )
    plt.xlim( right=150 )
    plt.xlabel( 'Lifespan', fontsize=16 )
    plt.ylabel( 'Normalized distribution', fontsize=16 )
    return fig


def plot_hist_overlay( ts, persfile=None, single=True, norm_it=True,
                       nbins=100, cell_type='new', vline=None ):
    """
    Plot the histogram figure for the RBC paper.

    ts -- List of arrays of times series of generator lifespans for
    each cell.

    persfile -- full path to a single persistence file (i.e. fsingle frame of single cell)

    single -- Cell to choose from list to compute histogram statistics
    on. Should be the same cell used in <persfile>. (Default=1,
    corresponds to new11 in ordered list (see below)).

    norm_it -- Toggle whether to return a normalized histogram.

    nbins -- number of bins.

    new -- New or Old cells.

    vline -- x-axis location of vertical dotted line. If None, no line is drawn.

    Note: Values used in RBC paper:

    new_hist_ts.pkl
    old_hist_ts.pkl
    
    new_110125-concatenated-ASCII_2000_1.txt
    old_120125-concatenated-ASCII_2000_1.txt
    """
    # for text object $\tau^*$ below
    from matplotlib.text import Text
    if cell_type == 'new':
        color = 'blue'
        ctype = 'new'
    elif cell_type == 'old':
        color = 'red'
        ctype = 'old'
    else:
        color = 'green'
        ctype = 'all'
    # compute stats for all cells
    out_all = plot_hist_all( ts, norm_it=norm_it )
    allx = out_all[0]
    ally = out_all[1]

    # compute stats for chosen single cell
    if single is not False:
        out = plot_hist_all( [ts[single]], norm_it=norm_it )
        nx = out[0]
        ny = out[1]
        pdf_ny = pdf( nx, ny )

    # now normalize everything by dividing by total area ( y --> PDF )
    pdf_ally = pdf( allx, ally )

    # output some stats 
    # print "\int { pdf_ally } = ", ((nx[1:]-nx[:-1]) * pdf_ally[:-1]).sum()
    # print "\int { pdf_ny } = ", ((nx[1:]-nx[:-1]) *pdf_ny[:-1]).sum()

    fig = plt.figure()
    ax = fig.gca()
    #    ax.set_xscale( 'log' )
    ax.set_yscale( 'log' )
    #ax.set_aspect( 1 )


    if single is not False:
        ax.plot( nx, pdf_ally, lw=3, c='g', label='Mean, all '+ctype+' cells' )
        ax.plot( nx, pdf_ny, lw=3, c='m', marker='^', ms=8,
                 label='Mean, single '+ctype+' cell' )
    if vline:
        ax.axvline( vline, linestyle=':', color='k' )
    # add a histogram for a single frame
    if persfile:
        ts = numpy.asarray( get_ts( persfile ), dtype=numpy.int )
        n, bins = numpy.histogram( ts, bins=nbins, range=(nx.min(),nx.max()) )
        ts_pdf = pdf( bins[:-1], n )

        #print 'ts_pdf', ((bins[1:] - bins[:-1])*ts_pdf).sum()
        width = bins[1]-bins[0]  #nx[1]-nx[0]
        ax.bar( bins[:-1], ts_pdf, width=width, color=color, alpha=0.5,
                label='Single frame distribution' )
        #ax.plot( bins[:-1], ts_pdf, marker='o', ms=6, lw=3, label='Single frame' )

    #ax.set_xticklabels( [str(int(x)) for x in xticks], fontsize=20 )
    plt.ylim( 0.00001, 0.1 )
    plt.xlim( right=150 )

    if vline:
        # add a \tau^* i nthe right spot
        tks, labels = plt.xticks()
        tks = list( tks )
        tks.append( vline )
        tks.sort()
        # find index of new vline tick
        loc = tks.index( vline )
        tau = Text( text='$\tau^{*}$' )
        new_labs = []
        # insert \tau into correct spot (Text( '\tau' ) doesn't seem to
        # work)
        for x in tks:
            if x == vline:
                if ctype == 'new':
                    L_text = r'$L_{new}$'
                elif ctype == 'old':
                    L_text = r'$L_{old}$'
                else:
                    L_text = r'$L_{all}$'
                new_labs.append( L_text )
            else:
                new_labs.append( str( int(x) ) )
        ax.set_xticks( tks )
        ax.set_xticklabels( new_labs, fontsize=12 )

    # now back to normal labeling and stuff
    plt.xlabel( 'Lifespan', fontsize=16 )
    plt.ylabel( 'Normalized distribution', fontsize=16 )

    plt.legend()
    #fig.show()
    return fig#, bins, ts_pdf  #, n, bins

def pdf( xi, yi ):
    """
    Normalize f(x) = y in terms of probability distribution
    functions. Thus, it should return f such that 

    \int_{min(xi)}^{max(xi)} {f(x_i)*(x_{i+1}-x_{i})   = 1
    """
    dx = xi[1:] - xi[:-1]
    integ = numpy.sum( yi[:-1] * dx )
    thepdf = yi/integ
    return thepdf


def interpolate( ts, bins ):
    yp = None
    xi = numpy.linspace( 0, bins[-1], 200 )
    yi = stineman_interp( xi, bins[:-1], ts, yp )
    # interpolate upper and lower error bars to get envelope. doesn't work with log scale.
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

