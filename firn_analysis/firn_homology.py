from pyTopTools import chomp_image as C
from pyTopTools import pyImage as P
import numpy as np
import matplotlib.pyplot as plt
import cPickle as pkl
import pp, time, os


#chomp_path = '/sciclone/data10/jberwald/CT_Firn_Samples/chomp_files2/'
#chomp_path = './_cubfile.cub'

def compute_stack_betti_numbers( ims, cubfile='/var/tmp/_cubfile.cub' ):
    """
    Compute the homology on a stack of images.
    """
    CI = P.ChompImage( ims )
    CI.get_cubical_corners()
    CI.cub2file( cubfile )
    CI.run_chomp( cubfile )
    CI.extract_betti()
    return CI

def firn_stacks( low, high, height, fpath, fprefix, cubpath,
                 step=None, cap=False, make_top=True, make_bottom=True ):
    """Build a sequence of image stacks, then compute the homology of the stack:

    ( start + height, start + step + height, start + 2*step+height, ...
    top - step + height )

    Returns a dictionary of betti numbers, keyed by stacks, key==(bottom, top).

    low, high :  min and max of the stack _bases_

    height : height of each stack (number of frames used to form the stack)

    Optional args:
    -------------

    step : steps to take in range. default=None will set step==height,
    so there are no overlapping blocks.

    """
    if step is None:
        step = height
    # stack images from bottom to top and record betti numbers
    betti = {}

    for base in range( low, high, step ):
        # create a sequence of images to stack
        ims = []
        # grab each image from base --> base+height
        for x in range( height ):
            image_num = base + x
            im = P.PyImage( fpath + fprefix + str( image_num ) + '.bmp' )
            im.bmp2array()
            ims.append( im.data )
        
        ims = np.array( ims, dtype=np.int )
        if cap:
            block = P.cap_block( ims, dtype=np.int, top=make_top, bottom=make_bottom)
        else:
            block = np.asarray( ims, dtype=np.int )
        # returns a ChompImage object
        print "  Computing betti numbers on the block ending at "+str(image_num)

        cubfile = cubpath + 'block_'+str(low)+'_'+str(height)
        if cap:
            if make_bottom:
                cubfile += '_bottom'
            if make_top:
                cubfile += '_top'
        cubfile += '.cub'
        b = compute_stack_betti_numbers( block, cubfile=cubfile )
        key = ( base, base+x )
        betti[ key ] = b.betti

    # save the betti numbers of each stack  
    savename = cubpath + fprefix + 'min' + str(low) +\
               '_max'+str(high) + '_block' + str(h)
    if cap:
        if make_top:
            savename += '_top'
        if make_bottom:
            savename += '_bottom'
    with open( savename + '.pkl', 'w' ) as fh:
        pkl.dump( betti, fh )    

    return betti


def make_betti_fig( stack_height, dim, fdir, prefix, suffix='',
                    fig=None, color='b', shape='o', size=6, label=None ):
    """
    Example:

    make_betti_figs( 40, 1, 'path/to/data/betti/' ),

    where the dirctory contains files of the form

    prefix + '_min*' + '_max*' + '_h*.pkl'

    suffix : default for uncapped blocks of height 40. change to 'top_bottom.pkl',
    'top.pkl' or 'bottom.pkl' for various caps.
    """
    # string version of height
    height_ = str( stack_height )
    if len( suffix ) > 0:
        suffix = height_ + '_' + suffix
    else:
        suffix = height_
    flist = os.listdir( fdir )
    # just keep the files for <height>
    betti_list = [ x for x in flist if x.endswith( suffix +'.pkl' ) ]
    
    nx = []
    # extract the base of each block
    for fname in betti_list:
        m = fname.find( 'min' )
        low = int( fname[ m+3:m+7 ] )
        nx.append( low )

    # assume for now that we're interested in the creation of tunnels
    # when we cap
    betti1 = []
    betti2 = []
    for b in betti_list:
        with open( fdir + b ) as fh:
            data = pkl.load( fh )
        # [0] needed since values returns a list regardless
        b1 = data.values()[0][1]
        b2 = data.values()[0][2]
        betti1.append( b1 )
        betti2.append( b2 )
    betti1 = np.array( betti1, dtype=np.int )
    betti2 = np.array( betti2, dtype=np.int )
    
    print "betti numbers " + suffix
    print "b1: ", betti1
    print "b2: ", betti2
    
    if dim == 1:
        bettis = betti1
    else:
        bettis = betti2

    if fig is None:
        fig = plt.figure(figsize=(8, 6)) 
    ax = fig.gca()
    ax.plot( bettis, color+shape+'-', ms=size, label=label )
    ax.set_xlabel( "Block number (height="+ height_ +")", fontsize=16 )
    ax.set_ylabel( r"# of generators", fontsize=16 )
    ax.set_xticklabels( nx )
    ax.set_xlim( -1, len(bettis) )
    #ax.set_ylim( min(bettis)-1, max(bettis)+1 )

    ax.legend()

    return fig, (nx,betti1,betti2)

def make_combined_betti_figs( nx, data1, data2, height, dim ):
    """Comparison of capped and uncapped blocks. 

    Plot the homology of the blocks on top, the differences below.

    nx : Base of each block for x-axis ticks

    data1, data2 : 1D numpy arrays of betti numbers across blocks.

         -- data1 : uncapped

         -- data2 : capped (if these are switched, the difference plot
    will be negative)
    """
    height_ = str( height )

    # create the current figure to plot with
    fig = plt.figure()

    # upper 2 rows, betti numbers
    ax1 = plt.subplot2grid( (3,2),(0,0), colspan=2,rowspan=2 )
    # bottom row, differences due to capping
    ax2 = plt.subplot2grid( (3,2),(2,0), colspan=2,rowspan=1 )

    ax1.plot( data1, 'bo-', ms=6, label='Uncapped' )
    ax1.plot( data2, 'gs-', ms=6, label='Capped' )
    #ax1.set_xlabel( "Block base (height="+ height_ +")", fontsize=16 )
    ax1.set_ylabel( r"$\beta_{"+str( dim )+"}$", fontsize=16 )
    ax1.set_xticklabels( [] )
    #ax.set_xlim( -1, len(bettis) )

    # differences
    assert len( data1 ) == len( data2 ), \
        "Number of blocks must be the same for data1 and data2"

    diffs = data1 - data2

    ax2.plot( diffs, 'rd-', ms=6, label='Uncapped - Capped' )
    ax2.set_xlabel( "Block base (height="+ height_ +")", fontsize=16 )
    ax2.set_ylabel( r"$\beta_{"+str( dim )+"}$", fontsize=16 )
    ax2.set_xticklabels( nx )
    ax2.set_xlim( -1, len(data1) )

    #plt.show()
    ax1.legend()
    ax2.legend()

    return fig
    

if __name__ == "__main__":


    ###################################
    #
    # Run experiments below
    #
    ###################################
    import sys

    start = time.time()

    #path = '/data/CT_Firn_Sample/output23-10-3/'
    path = '/sciclone/data10/jberwald/CT_Firn_Samples/output23-10-3/'

    #chomp_path = '/data/CT_Firn_Sample/debug/'

    prefix = 'K09b-23-10-'

    #===========================
    # Build stacks or stacks with caps
    #===========================
    if 0:

        # change these two values
        #caps = True
        #caps = False  # just run both in series
        h = int( sys.argv[2] )

        # input from shell script
        low = int( sys.argv[1] )
        high = low + h

        print "height = ", h
        print low, high

        # do the homology computations

        print "Building blocks..."
        betti_blocks = '/sciclone/data10/jberwald/CT_Firn_Samples/firn_blocks_June2013/'
        betti_caps = '/sciclone/data10/jberwald/CT_Firn_Samples/firn_caps_June2013/'

        #Bblock = firn_stacks( low, high, h, path, prefix, betti_blocks )
        #Bcap = firn_stacks( low, high, h, path, prefix, betti_caps, cap=True ) 
        Bcap_top = firn_stacks( low, high, h, path, prefix, betti_caps, cap=True, make_bottom=False )  
        Bcap_bottom = firn_stacks( low, high, h, path, prefix, betti_caps, cap=True, make_top=False )  
       
        print "Total time: ", time.time() - start

    if 1:

        betti_blocks = '/data/CT_Firn_Sample/firn_blocks_June2013/'
        betti_caps = '/data/CT_Firn_Sample/firn_caps_June2013/'
        #dims = [0,1,2] #[0,1,2]
        stack_height = [60]

        dim = 1
        for h in stack_height:
            print "stack height", h
            fig, betti = make_betti_fig( h, dim, betti_blocks, prefix, label=r'$\beta_1$, uncapped' )
            fig, btop = make_betti_fig( h, dim+1, betti_caps, prefix, fig=fig, color='g', shape='s',
                                        suffix='top', label=r'$\beta_2$, top' )
            fig, bbottom = make_betti_fig( h, dim+1, betti_caps, prefix, fig=fig, color='r', shape='s',
                                           suffix='bottom', label=r'$\beta_2$, bottom' )

            fig, btb = make_betti_fig( h, dim+1, betti_caps, prefix, fig=fig, shape='d', color='k',
                                        suffix='top_bottom', label=r'$\beta_2$, top & bottom' )
            # fig, b2 = make_betti_fig( h, dim+1, betti_caps, prefix, fig=fig, color='r', 
            #                           suffix=suffix, shape='s', label=r'$\beta_2$, capped' )
            
            # fig.savefig( '/data/CT_Firn_Sample/figures/blocks_caps/'+\
            #              'betti_h'+str(h)+'_d'+str(dim)+'.pdf' )
            fig.show()
            

            # compute  \beta_2^{top/bottom} - ( \beta_2^{top}  +\beta_2^{bottom} -2\beta_2 )
            # external tunnels due to one-sided capping
            bx = (btop[2] - betti[2]) + (bbottom[2] - betti[2])
            # full tunnels = \beta_2^{top/bottom} - bx - \beta_2
            b2full = btb[2] - bx - betti[2]

            
            # now, compute tunnels lost to capping
            nx = betti[0] # x-axis

            # full tunnels == \beta_2^{top/bottom} - ( \beta_2^{top}  +\beta_2^{bottom} )
            # d2 = b2[1]
            # fig_comb = make_combined_betti_figs( nx, d1, d2, h, dim )
            # fig_comb.show()
