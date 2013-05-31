import numpy as np
import matplotlib.pyplot as plt


def persdia_sub_super( filename, **kwargs ):
    """
    filename -- path to .npy or .txt file containing persistence intervals

    birth and death info in column form
    
    """
    args = { 'ms' : 2 }
    args.update( kwargs )

    
    # set a flag in case the only death is infinite; default is false
    maxd_is_inf = False 
    normal_gens = False

    # extract birth and death indices
    try:
        ints = np.load( filename )
    except IOError:
        ints = np.loadtxt( filename )
    births = ints[:,0]
    deaths = ints[:,1]

    # extract indices of those intervals which persist throughout
    sub_ints = np.where( deaths == -1 )[0]
    super_ints = np.where( deaths == -2 )[0]

    # extract indices of those intervals which die
    try:
        n1 = np.where( deaths != -1 )[0]
        n2 = np.where( deaths != -2 )[0]
        normal_ints = np.intersect1d( n1, n2 )
        normal_gens = True  # we have normal (non-infinite gens)
    except:
        normal_ints = []
        print 'Only infinite generators on this homology level.' 

    # extract maximum death time encountered and min birth times
    #minb = births.min()  #min( births ) % should always be >= 1
    maxbirth = births.max() #max( births )
    maxd = deaths.max() # max( deaths )
    if maxd < 0:
        maxd_is_inf = True;
        maxd = births.max() #max( births )

    # start drawing the figure
    fig = plt.figure()
    ax = fig.gca()
    
    # we always plot these
    ax.plot( births[ normal_ints ], deaths[ normal_ints ], 'bo', **args  );

    # create diagonal
    diag = [0, maxd+2]
    ax.plot(diag, diag, 'g-')

    # set super-level infinite gens to min - 1 (so if min val in
    # matrix is 2, set inf to 1 )
    max_axis = max( [ maxd, maxbirth ] )
    infsub_vec = (max_axis + 1) * np.ones_like( sub_ints )
    #%infsuper_vec = (minb -1 ) * ones(size(super_ints));
    infsuper_vec = np.zeros_like( super_ints )

    # create diag, from equally spaced points (20 steps)
    diag = np.linspace( 0, (1.1)*(max_axis+2), (1.1)*(max_axis+2)/20 )
    npts = len( diag );

    # plot diag and horizontal lines
    ax.plot(diag,diag,'g-');
    ax.plot( npts * [0], diag, 'k--' )
    ax.plot( diag, npts * [0], 'k--' )

    # plot infinite stuff
    ax.plot( births[ sub_ints ], infsub_vec , 'rd', lw=1, **args )
    ax.plot( births[ super_ints ], infsuper_vec , 'rd', **args )

    # axis lines, dashed
    ax.hlines( 0, 0, max_axis+1, linestyle='dashed' )
    ax.vlines( 0, 0, max_axis+1, linestyle='dashed' )

    # set axis limits
    ax.set_xlim( [ -0.1*maxd, (1.1)*(max_axis+1) ] )
    ax.set_ylim( [ -0.1*maxd, (1.1)*(max_axis+1) ] )

    ax.set_title(filename)
    ax.set_xlabel('birth')
    ax.set_ylabel('death')

    return fig
