import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import cPickle as pkl

# local directory
import gaussian as gauss
import rbc_npy2Perseus as rp
import rbc_perseus as pers
import rbc_postprocess as rpost
import rbc_histogram as rh
import robust_generators_hist as rg

# SET GLOBAL 'SHOW PLOT' VARIABLE
global show_plts
show_plts = False


def plot_frame( frame, outname=None ):
    """
    Use pylab.imshow() to plot a frame (npy array). 
    Note: The boundary should be removed, thus there will not be a halo

    Default frame is new_110125-concatenated-ASCII_1050.npy (see
    make_figs function).
    """
    data = np.load( frame )
    fig = plt.figure()
    #plt.title("RBC frame")
    plt.imshow(data)
    plt.colorbar()
    if show_plts:
        plt.show()
    if outname:
        fig.savefig( outname )
    return fig


################################
# GAUSSIAN FIGS
################################

def plot_gauss_bump( noise=None, clip=True ):
    """
    Returns an image of a matrix with two Gaussian bumps combined (the
    small and large bump examples).
    """
    G = gauss.gauss_bump( noise=noise )
    if clip:
        gm = np.ma.masked_less_equal( G, 1 )
        
    fig = plt.figure( frameon=False )
    plt.axes( frameon=False )
    ax = fig.gca()
    ax.set_xticks([])
    ax.set_xticklabels([])
    ax.set_yticks([])
    ax.set_yticklabels([])

    cax = ax.imshow( gm )
    fig.colorbar( cax )
    if show_plts:
        fig.show()

    # return the figure and the data (for sublevel use)
    return fig, gm

def plot_sublevel( G, level ):
    """
    Levels 2, 10, 18 used in paper.
    """
    fig, S = gauss.sublevel( G, level )
    return fig, S

def plot_gauss_sublevels( level ):
    # no noise
    data = gauss.gauss_bump()
    fig, G = gauss.sublevel( data, level, show_plt=show_plts )

def create_persfile( smooth_data,
                     noise=None,
                     scale=1000 ):
    """
    smooth_data : non-noisy Gaussian. Used to determine proper outline
    for clipping noisy Gaussian.

    noise : noise level epsilon (set to 0.05 in paper)

    scale : for proper resolution in Perseus, which only accept
    integers for sublevel set computations.

    Saves the perseus-readable file (sparse cubical format) to disk (see
    below).
    """
    # grab the original file
    B = smooth_data

    if noise:
        A = gauss.gauss_bump( noise=noise )
        # find where to clip the noisy surface 
        gauss.clip_below( B, 1 ) # in-place, set elements to zero
        w = np.where( B == 0 )
        # now both smooth and noisy surfaces are clipped to zero outside
        # the same boundary.
        A[w] = 0
        
        # scale for additional resolution
        A *= scale

        # might as well save the file for posterity
        noise_level = str( noise )
        # remove the decimal (split on '.', then join the list with '.' removed)
        noise_level = ''.join( noise_level.split( '.' ) )
        sname = './data_for_figs/gauss_Gnoise'+noise_level
    else:
        sname = './data_for_figs/gauss_smooth'
        A = B

    # save the noisy matrix
    print "Saving ", sname
    np.save( sname+'.npy', A )
    # read it and output the proper Perseus format
    persname = sname + '_pers'
    rp.write_sparse_file( sname+'.npy',
                          sname+'_pers' )
    return persname

################################
# END GAUSSIAN FIGS
################################


    
################################
# RUN PERSEUS
################################
def run_perseus( fname, outname ):
    pers.perseus( fname, outname )


################################
# SUBLEVEL FIGS
################################

def plot_frame_mask_zero( frame, nx=203, ny=198 ):
    """
    Use pylab.imshow() to plot a frame. Mask the elements of the image
    matrix that are zero.

    nx and ny are typical dimensions of the cell frames. Some cells
    are different so these will have to be adjusted accordingly.
    """
    # This allows frame to be any data, such as a symmetric 2D
    # Gaussian. 
    try:
        frame.resize((nx,ny))
    except:
        pass

    cdict = {'red': ((0., 1, 1),
                     (0.05, 1, 1),
                     (0.11, 0, 0),
                     (0.66, 1, 1),
                     (0.89, 1, 1),
                     (1, 0.5, 0.5)),
            'green': ((0., 1, 1),
                      (0.05, 1, 1),
                      (0.11, 0, 0),
                      (0.375, 1, 1),
                      (0.64, 1, 1),
                      (0.91, 0, 0),
                      (1, 0, 0)),
            'blue': ((0., 1, 1),
                     (0.05, 1, 1),
                     (0.11, 1, 1),
                     (0.34, 1, 1),
                     (0.65, 0, 0),
                     (1, 0, 0))}

    my_cmap = colors.LinearSegmentedColormap('my_colormap',cdict,256)
    fig = plt.figure( frameon=False, dpi=160 )
    # for transparency
    fig.patch.set_alpha( 0.0 )
    ax = fig.gca()
    ax.set_xticks( [] )
    ax.set_yticks( [] )
    ax.set_axis_off() # turn off the axes 
    im = ax.imshow( frame, cmap=my_cmap )
    cbar = plt.colorbar( im, shrink=0.6 )
    fig.show()
    return fig
    
def plot_sublevel_set( frame, height, bndfile=None, persfile=None,
                       nx=203, ny=198, save=False,
                       transparent=True, thresh=None ):
    """
    Plot sublevel set for an array representing an intensity function.

    frame -- path to array on disk or Numpy array.

    height -- function height at which to take sublevel set.

    Returns figure object
    """
    h = height
    # text file
    try:
        try:
            data = np.loadtxt( frame )
            # np file
        except ValueError:
            data = np.load( frame )
    except:
        # This is for a general array that is already rectangular
        # and needs no help.
        data = frame
        nx,ny = data.shape

    print "data", data.shape
    # if we read in a data file we might need a boundary file to go with it.
    if bndfile:
        bnd = np.loadtxt( bndfile )
        print "boundary", bnd.shape
        try:
            data = bnd.ravel() * data
        except ValueError:
            data = bnd * data
        data.resize( (nx,ny) )

    # make an array to hold (R,G,B,A) data at each pixel
    G = np.zeros((nx,ny,4), dtype=int)

    # if not thresh:
    #     thresh = 0.9
    
    #temp = data.copy()
    G[ np.where( data > int(h) ) ] = [1,1,1,0]
    # the sublevel set
    G[ np.where( data <= int(h) ) ] = [0,0,160,1]
    # everything outside 
    G[ np.where( data == 0.0 ) ] = [1,1,1,0]
    # outside majority of Gaussian peak
    if thresh:
        G[ np.where( data <= thresh ) ] = [1,1,1,0]
   
    # now plot stuff
    fig = plt.figure( frameon=False, dpi=160 )

    # make things transparent?
    if transparent:
        fig.patch.set_alpha( 0.0 )

    ax = fig.gca()
    ax.set_frame_on( False )
    #ax.set_title('RBC ' + output)
    # PLOT THE MATRIX OF VALUES
    print "plotting sublevel set..."
    ax.imshow( G )
    ax.set_xticks( [] )
    ax.set_yticks( [] )
    if show_plts:
        fig.show()
    return fig

def sublevel_sequence():
    heights = [ 900,1050,1110 ]
    frame = 'new11_frame2000.txt'
    bnd = 'boundary_Nov_new110125'
    figs = []
    for h in heights:
        figs.append( plot_sublevel_set( frame, h, bndfile=bnd, mask=True ) )
    return figs

################################
# END SUBLEVEL FIGS
################################

################################
# HISTOGRAMS
################################
            
def plot_hist_figure( ts, persfile=None, single=1, norm_it=True,
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
        ts = np.asarray( get_ts( persfile ), dtype=np.int )
        n, bins = np.histogram( ts, bins=nbins, range=(nx.min(),nx.max()) )
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
    fig.show()
    return fig, bins, ts_pdf  #, n, bins

def plot_hist( fname, nbins=None, scale=1.0, color='blue',
               gaussian=False, sigma=1.0, normed=False, fontsize=14 ):
    """
    Plot a histogram of generator lifespan along the diagonal.

    fname -- full path to perseus output file.
    """
    # for normpdf() function
    import matplotlib.mlab as mlab
    
    ts = get_ts ( fname )
    # the (almost) infinite generator overwhelms the plot
    ts = ts[:-1]

    print ts
    if scale:
        ts = numpy.asarray( ts, dtype=numpy.float )
        ts /= scale

    # plot the histogram
    fig = plt.figure()
    ax = fig.gca()
    # the histogram of the data
    if not nbins:
        nbins = ts.max()-ts.min()
    n, bins, patches = ax.hist( ts, bins=nbins, normed=normed,
                                facecolor=color, alpha=0.75)

    if gaussian:
        mu = 0
        bincenters = 0.5 * ( bins[1:] + bins[:-1] )
        # add a 'best fit' line for the normal PDF
        y = mlab.normpdf( bincenters, mu, sigma)
        l = ax.plot( bincenters, y, 'r--', linewidth=2 )

    # xticks = [ int( tk ) for tk in ax.get_xticks() ]
    # yticks = [ int( tk ) for tk in ax.get_yticks() ]
    # ax.set_xticklabels( xticks, fontsize=fontsize )
    # ax.set_yticklabels( yticks, fontsize=fontsize )

    #ax.set_title( r'Distribution of generator lifespans along diagonal' )
    ax.set_xlabel( r"Lifespan", fontsize=fontsize )
    ax.set_ylabel( r"Number of generators ($\beta_1$)", fontsize=fontsize )
    #ax.grid( True )

    ax.set_xlim( 0,21 )
    
    plt.show()
    return fig, ts, bins

def midrange_stats( fname=None, new=True, stats=False, ngens=False ):
    """
    For the midrange generators in each frame in each cell, find a
    tuple = ( median, mean, std ). 
    """
    newlist = ['new_10' , 'new_110125', 'new_130125', 'new_140125', 'new_3',
               'new_4', 'new_40125', 'new_50125', 'new_6', 'new_60125', 'new_9']
    oldlist = ['old_100125', 'old_120125', 'old_15', 'old_2', 'old_4000', 'old_4001',
               'old_5',  'old_50125',  'old_6',  'old_7',  'old_8',  'old_9',  'old_90125']

    if fname:
        print "\nreading " + fname + " ... "
        if new:
            cell_list = newlist
        else:
            cell_list = oldlist
        with open( fname ) as fh:
            gens = pkl.load( fh )
    else:
        # previous versions, static file names
        if new:
            print "computing New means...\n"
            cell_list = newlist
            with open( '/Users/jberwald/Dropbox/Projects/rbc/pyRBC/data/new_midrange_means_normed_eps30.pkl' ) as fh:
                gens = pkl.load( fh )
        else:
            print "computing Old means...\n"
            cell_list = oldlist
            with open( '/Users/jberwald/Dropbox/Projects/rbc/pyRBC/data/old_midrange_means_normed_eps30.pkl' ) as fh:
                gens = pkl.load( fh )
        
    # store (median, means, std) in list, within dictionary keyed by cell name
    data_dict = {}

    if stats:
        # compute the stats
        for i, cell in enumerate( cell_list ):
            ci = gens[i]
            # loop over each frame's (birth,death) coords
            # store stats for lifespan for each in a list
            stats = []
            for x in ci:
                try:
                    diff = x[:,1] - x[:,0]
                    data = ( numpy.median( diff ), diff.mean(), diff.std() )
                except TypeError:
                    data = ( 0., 0., 0. )
                stats.append( data )
            data_dict[ cell ] = numpy.asarray( stats )
    elif ngens:
        # concatenate the data and return
        for i, cell in enumerate( cell_list ):
            ci = gens[i]
            numgens = []
            for x in ci:
                try:
                    numgens.append( len( x ) )
                except TypeError:
                    numgens.append( 0 )
            data_dict[ cell ] = numpy.array( numgens, dtype=int )

    return data_dict

def find_maxes():
    # find the old maxes. do these in sequence to avoid killing
    # the disk with thousands of minute searches (sigh).
    old_prefix = '/data/PerseusData/PerseusOutput/original/2d_sparse/Old/'
    old_dirs= os.listdir( old_prefix )
    old_cells = [ old_prefix + c + slash for c in old_dirs ]
    print "finding maxima for old cells..."
    with timer.Timer():
        old_maxes = all_maxes( old_cells )
        with open( './data/old_maxes.pkl', 'w' ) as fh:
            pkl.dump( old_maxes, fh )

    # now find the new maxes
    new_prefix = '/data/PerseusData/PerseusOutput/original/2d_sparse/New/'
    new_dirs= os.listdir( new_prefix )
    new_cells = [ new_prefix + c + slash for c in new_dirs ]
    print "finding maxima for new cells..."
    with timer.Timer():
        new_maxes = all_maxes( new_cells )
        with open( './data/new_maxes.pkl', 'w' ) as fh:
            pkl.dump( new_maxes, fh )


################################
# HISTOGRAMS
################################
            
################################
# END PERSISTENCE DIAGRAMS
################################
            
def midrange_timeseries( old_prefix, new_prefix, eps1, eps2, normed=False, plot_fig=False ):
    """
    Time series of midrange (==robust) generators.
    """
    # if normed:
    #     old_prefix = '/Users/jberwald/Dropbox/Projects/rbc/pyRBC/data/old_midrange_means_normed_eps'
    #     new_prefix = '/Users/jberwald/Dropbox/Projects/rbc/pyRBC/data/new_midrange_means_normed_eps'
    # else:
    #     old_prefix = '/Users/jberwald/Dropbox/Projects/rbc/pyRBC/data/old_midrange_eps'
    #     new_prefix = '/Users/jberwald/Dropbox/Projects/rbc/pyRBC/data/new_midrange_eps'

    
    old_fname = old_prefix + str(eps1) + "_" + str(eps2) + ".pkl"
    new_fname = new_prefix + str(eps1) + "_" + str(eps2) + ".pkl"
    # dictionary keyed=cells, values=array of number of midrange
    # generators for each frame of cell
    old_mr = midrange_stats( fname=old_fname, ngens=True )
    new_mr = midrange_stats( fname=new_fname, ngens=True )

    if plot_fig:
        figs = []
        oldkeys = old_mr.keys()
        newkeys = new_mr.keys()
        # choose a couple of cells to plot
        for i in [7,8]: #range(3):
            oldname = oldkeys[i]
            newname = newkeys[i]
            oldvals = old_mr[ oldkeys[i] ]
            newvals = new_mr[ newkeys[i] ]
            fig = plot_midrange_ts( newvals, oldvals, skip=10, means=True,
                                    fontsize=14 )
            ax = fig.gca()
            #title = "Old cell vs New cell ("+str(i)+")"
            #ax.set_title( title )

            # adjust to tighten the bbox and make sure the xlabel fits
            # on the figure
            fig.subplots_adjust( bottom=0.14, left=0.07, right=0.96, top=0.95 )
                
            fig.show()
            # fig.savefig( '/Users/jberwald/Dropbox/Projects/rbc/pyRBC/data/'+\
            #              'midrange_ts_'+oldname+'-'+newname+'eps'+\
            #              str(eps1)+'_'+str(eps2)+'.png',
            #              dpi=400, transparent=True )    
            
        
            # now add means to the plot
            #new_mean = newvals.mean()
            #old_mean = oldvals.mean()

            # ax = fig.gca()
            # ax.axhline( new_mean, color='b', linestyle='--', lw=2, label="New cell mean" )
            # ax.axhline( old_mean, color='r', linestyle='--', lw=2, label="Old cell mean" )
            # ax.legend()
            figs.append( fig )
        return figs
    else:
        return old_mr, new_mr


def make_all_figs():
    """
    Comment/uncomment functions according to desired figs.

    The data for the figs is stored within the figs folder.

    The fig names are aligned with the current version of the paper
    (June 15, 2013).

    One must create a figs/ directory at the local level.
    """
    # Change value up here if desired
    # noise level for gaussian bump example
    height = 1050  # for sublevel set
    noise = eps = 0.05 # for Gaussian
    alpha = 0.4  # transparency level

    data_prefix = 'data_for_figs/'
    save_prefix = './figs/'
    # store and return all fig objects for inspection
    figs = {}

    #======================
    # FIGURE 1
    #======================
    if 1:
        print "Making Figure 1..."
        # single frame of one cell
        frame_path = data_prefix + 'new11_frame2000.npy'
        fig1 = plot_frame( frame_path )
        fig1.savefig( save_prefix + 'new11_frame2000_cbar.pdf' )

        # sublevel set
        fig1b = plot_sublevel_set( frame_path, height )
        fig1b.savefig( save_prefix + 'new11_frame2000_sub1050.pdf' )
        figs['cell_sublevel'] = fig1b

        #======================
        # FIGURE 2 & 3
        #======================
        print ""
        print "Making Figure 2 & 3..."
        print "  Gaussians and sublevels (note, must use Sage for 3d plots)..."
        # smooth gaussian (2d)
        fig2, gsmooth = plot_gauss_bump()
        fig2.savefig( save_prefix + 'gauss_peak_cbar.pdf' )
        figs['gauss2d_smooth'] = fig2

        # noisy gaussian (2d)
        fig3, gnoise = plot_gauss_bump( noise=noise )
        fig3.savefig( save_prefix + 'gauss_peak_cbar.pdf' )
        figs['gauss2d_noise'] = fig3

        # sublevels for gaussians above these require post-processing in
        # order to create the sequences Figures 2(b) and 3(b). Performed
        # using Inkscape/Gimp
        fig4, gs2 = plot_sublevel( gsmooth, 2 )
        fig4.savefig( save_prefix + 'gsmooth_sublevel2.pdf' )
        fig5, gs10 = plot_sublevel( gsmooth, 10 )
        fig5.savefig( save_prefix + 'gsmooth_sublevel10.pdf' )
        fig6, gs18 = plot_sublevel( gsmooth, 18 )
        fig6.savefig( save_prefix + 'gsmooth_sublevel18.pdf' )
        figs['gs2'] = fig4
        figs['gs10'] = fig5
        figs['gs18'] = fig6

        fig7, gn2 = plot_sublevel( gnoise, 2 )
        fig7.savefig( save_prefix + 'gnoise_sublevel2.pdf' )
        fig8, gn10 = plot_sublevel( gnoise, 10 )
        fig8.savefig( save_prefix + 'gnoise_sublevel10.pdf' )
        fig9, gn19 = plot_sublevel( gnoise, 19 )
        fig9.savefig( save_prefix + 'gnoise_sublevel19.pdf' )
        figs['gn2'] = fig7
        figs['gn10'] = fig8
        figs['gn19'] = fig9


        ##########################
        # SMOOTH GAUSSIAN DIAGRAM
        ##########################
        print "  Persistence diagrams..."
        # Perseus diagram for smooth Gaussian
        gdata = gauss.gauss_bump()
        persname = create_persfile( gdata, scale=1 )
        # strip '_pers' of the end
        outname = persname[:-5]
        run_perseus( persname+'.txt', outname )

        # Plot the smooth Gaussian persistence diagram
        fig10 = rpost.plot_diagram_std( outname + '_1.txt', scale=1, 
                                        show_fig=False, color='r', shape='s' )
        if show_plts:
            fig10.show()
        figname = save_prefix + 'gauss_smooth_dia.pdf'
        # figname =  outname[:-6] + 'dia.png'
        # s = figname.split( '/' )
        # idx = s.index( 'data_for_figs' )
        # s[ idx ] = 'figs'
        # figname = '/'.join( s )
        fig10.savefig( figname, dpi=200 )
        figs['gspers'] = fig10

        # The Perseus computation is not instantaneous for the noisy
        # Gaussian. Set to 'if 1' to run Perseus on the noisy version.
        if 1:
            ##########################
            # NOISY GAUSSIAN DIAGRAM
            ##########################
            print ""
            print "\tComputing persistence diagram for noisy Gaussian. This takes a little bit..."
            # Create persistence diagrams for noise level == \epsilon=0.05
            # NOTE: This plots to overlayed diagrams in figure 3
            eps = 0.05
            scale = 1000
            persname = create_persfile( gsmooth, noise=eps, scale=scale )
            # strip '_pers' of the end
            outname = persname[:-5]
            run_perseus( persname+'.txt', outname )

            # # Plot the noisy Gaussian persistence diagram
            # fig10b = rpost.plot_diagram_std( outname + '_1.txt', scale=,
            #                                  show_fig=False, color='r', shape='s' )
            # Uncomment and comment above fig10b if plotting noisy diagram without overlay
            # fig10 = None  
            fig11 = rpost.plot_diagram_std( outname + '_1.txt', scale=scale,
                                            show_fig=False, fig=fig10 )
            if show_plts:
                fig11.show()
            figname = save_prefix + 'gauss_noise'+str(eps)+'_dia.pdf'
            # figname = persname + '_dia.png'
            # s = figname.split( '/' )
            # idx = s.index( 'data_for_figs' )
            # s[ idx ] = 'figs'
            # figname = '/'.join( s )
            fig11.savefig( figname, dpi=200 )
            figs['gnpers'] = fig11

        #======================
        # FIGURE 4
        #======================
        print ""
        print "Making Figure 4..."
        # OLD CELL
        #==================
        frame_path = data_prefix + 'old12_frame2000.npy'
        fig12= plot_frame( frame_path )
        fig12.savefig( save_prefix + 'old12_frame2000.pdf' )
        figs['cell_old'] = fig12

        # Sublevel set
        subfigs = []
        for h in [1050,1450]:
            fig12b = plot_sublevel_set( frame_path, h )
            fig12b.savefig( save_prefix + 'old12_frame2000_sub'+str(h)+'.pdf' )
            subfigs.append( fig12b )
        figs['cell_old_sub'] = subfigs

        # PD in Figure 4
        persname =  frame_path[:-4]+'_pers'
        outname = frame_path[:-4]
        rp.write_sparse_file( frame_path,
                              frame_path[:-4]+'_pers' )
        run_perseus( persname + '.txt', outname )
        fig12c = rpost.plot_diagram_std( outname + '_1.txt', scale=1, 
                                        show_fig=False, color='b', shape='o' )
        fig12c.savefig( save_prefix + 'old12_frame2000_dia.pdf' )
        figs['cell_old_pd'] = fig12c

        # PD with regions (OLD)
        #========================
        print ""
        print "  Making OLD cell persistence diagram with regions..."
        lines = [1050,1450]
        fig12d = rpost.plot_diagram_regions( outname + '_1.txt', lines=lines, zoom=False )
        # fill in regions with alpha color
        ax = fig12d.gca()
        for line in lines:
            rect = matplotlib.patches.Rectangle((0,line),line,
                                                ax.get_ylim()[1] - line,
                                                color='#00C5CD', alpha=0.4)
            ax.add_patch( rect )
        ax.set_xlim( (lines[0]-200, lines[1]+200) )
        ax.set_ylim( (lines[0]-200, lines[1]+200) )
        ax.set_autoscale_on( False )

        fig12d.savefig( save_prefix + 'old12_frame2000_dia_regions.pdf' )
        figs['cell_old_regions'] = fig12d


        # NEW CELL
        #===================
        # Figure 1 has already produced the heat map for new11_frame2000
        # Sublevel set
        frame_path = data_prefix + 'new11_frame2000.npy'
        subfigs = []
        for h in [1050,1110]:
            fig13 = plot_sublevel_set( frame_path, h )
            fig13.savefig( save_prefix + 'new11_frame2000_sub'+str(h)+'.pdf' )
            subfigs.append( fig13 )
        figs['cell_old_sub'] = subfigs

        # PD in Figure 4
        persname =  frame_path[:-4]+'_pers'
        outname = frame_path[:-4]
        rp.write_sparse_file( frame_path,
                              frame_path[:-4]+'_pers' )
        run_perseus( persname + '.txt', outname )
        fig13c = rpost.plot_diagram_std( outname + '_1.txt', scale=1, 
                                         show_fig=False, color='b', shape='o' )
        fig13c.savefig( save_prefix + 'new11_frame2000_dia.pdf' )
        figs['cell_new_pd'] = fig13c



        # PD with regions (NEW)
        #=========================
        print ""
        print "  Making NEW cell persistence diagram with regions..."
        lines = [1050,1110]
        fig13d = rpost.plot_diagram_regions( outname + '_1.txt', lines=lines, zoom=False )
        # fill in regions with alpha color
        ax = fig13d.gca()
        for line in lines:
            rect = matplotlib.patches.Rectangle((0,line),line,
                                                ax.get_ylim()[1] - line,
                                                color='#00C5CD', alpha=0.4)
            ax.add_patch( rect )
        ax.set_xlim( (lines[0]-200, lines[1]+200) )
        ax.set_ylim( (lines[0]-200, lines[1]+200) )
        ax.set_autoscale_on( False )
        fig13d.savefig( save_prefix + 'new11_frame2000_dia_regions.pdf' )
        figs['cell_new_regions'] = fig13d
                                         
    #=====================
    # FIGURE 5
    #=====================

    # All-lifespan histograms
    #==========================
    print ""
    print "Making Figure 5..."
    old = data_prefix + 'old_hist_ts.pkl'
    new = data_prefix + 'new_hist_ts.pkl'
    #  new = fdir + 'new_hist_ts.pkl'
    #old = fdir + 'old_hist_ts.pkl'

    print "  Reading timeseries of NEW lifespans..."
    with open(new) as fh:
        A = pkl.load(fh)
        
    # print "\tMaking histogram..."
    # fig14a = rh.plot_hist( A, nbins=200, alpha=0.7 )
    # fig14a.savefig( save_prefix + 'new_hist_all.pdf' )

    print "  Reading timeseries of OLD lifespans..."
    with open(old) as fh:
        B = pkl.load(fh)

    # print "\tMaking histogram..."
    # fig14b = rh.plot_hist( B, nbins=200, color='r', alpha=0.7 )
    # fig14b.savefig( save_prefix + 'old_hist_all.pdf' )

    # Zoomed figures with single cell and overlay
    #=============================================
    print ""
    print "Making zoomed, overlayed histograms..."
    print "\tSometimes a 'divide-by-zero' error occurs. Don't panic, it's just due to "\
        "the interpolation algorithm not finding generators in certain histogram bins."
    print ""
    new_cell = data_prefix + 'new_110125-concatenated-ASCII_2000_1.txt'
    old_cell = data_prefix + 'old_120125-concatenated-ASCII_2000_1.txt'
    fig15a = rh.plot_hist_overlay( A, persfile=new_cell, nbins=500, vline=70 )
    fig15b = rh.plot_hist_overlay( B, persfile=old_cell, nbins=500, vline=55, cell_type='old' )
    fig15a.savefig( save_prefix + 'new_hist_overlay.pdf' )
    fig15b.savefig( save_prefix + 'old_hist_overlay.pdf' )

    #==========================
    # FIGURE 6
    #==========================
    print ""
    print "Making Figure 6..."
    bt_figs = rg.birth_time_hist( data_prefix, save_prefix )
    figs['birth_times'] = bt_figs
        
    return figs

if __name__ == "__main__":

    # import argparse
    # import sys
    
    # ##############################
    # #  PARSE COMMAND LINE ARGS
    # ##############################
    # parser = argparse.ArgumentParser(description="Module for producing figures in RBC paper. "\
    #                                  "Process input for rbc_paper_figs.py." )
    # parser.add_argument( "--cell", help="Produce Figure 1: heat map of single frame of a cell. "\
    #                      "See plot_frame() for details", action="store_true" )
    # parser.add_argument( "--gauss_bump", help="Plot the heat map for the two-bump Gaussian in Figure 2: "\
    #                      "See plot_gauss_bump() for details", action="store_true" )
    # parser.add_argument( "--gauss_bump_noise", help="Plot the heat map for the noisy two-bump Gaussian in Figure 3:"\
    #                      "See plot_gauss_bump() for details", action="store_true" )
    # parser.add_argument( "-m", "--matlab", help="Input data is in matlab matrices.",
    #                      action="store_true" )
    
    # args = parser.parse_args()

    # if len( sys.argv ) == 1:
    #     parser.print_help()
    # else:


    # plot some figures!
    F = make_all_figs()
        

