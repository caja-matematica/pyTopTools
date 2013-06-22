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
    block = P.build_block( ims )

    # cubfile = chomp_path + prefix[:-1] + \
    #     '_b' + str( base ) + \
    #     '_h' + str( height ) + '.cub'

    CI = P.ChompImage( ims )
    CI.get_cubical_corners()
    CI.cub2file( cubfile )
    CI.run_chomp( cubfile )
    CI.extract_betti()
    return CI

def firn_stacks( bottom, top, height, fpath, fprefix, cubpath, step=None, cap=False ):
    """
    Build a sequence of image stacks, then compute the homology of the stack:

    ( start + height, start + step + height, start + 2*step+height, ...
    top - step + height )

    Returns a dictionary of betti numbers, keyed by stacks, key==(bottom, top).

    Optional args:
    -------------

    step : steps to take in range. default=None will set step==height,
    so there are no overlapping blocks.
    """
    if step is None:
        step = height
    # stack images from bottom to top and record betti numbers
    betti = {}

    for base in range( bottom, top, step ):
        # create a sequence of images to stack
        ims = []
        # grab each image from base --> base+height
        for x in range( height ):
            image_num = base + x
            im = P.PyImage( fpath + fprefix + str( image_num ) + '.bmp' )
            im.bmp2array()
            ims.append( im.data )

        if cap:
            block = P.cap_block( ims )
        else:
            block = P.build_block( ims )
        # returns a ChompImage object
        print "  Computing betti numbers on the block ending at "+str(image_num)

        cubfile = cubpath + 'block_'+str(bottom)+'_'+str(height)+'.cub'
        b = compute_stack_betti_numbers( block, cubfile=cubfile )
        key = ( base, base+x )
        betti[ key ] = b.betti

    return betti


def make_betti_fig( stack_height, dim, fdir, prefix,
                    fig=None, color='b', shape='o', size=6, label=None ):
    """
    Example:

    make_betti_figs( 40, 1, 'path/to/data/betti/' ),

    where the dirctory contains files of the form

    prefix + '_min*' + '_max*' + '_h*.pkl'
    """
    # string version of height
    height_ = str( stack_height )
    flist = os.listdir( fdir )
    # just keep the files for <height>
    betti_list = [ x for x in flist if x.endswith( height_ +'.pkl' ) ]

    print "Number of blocks", len(betti_list)
    
    nx = []
    # extract the base of each block
    for fname in betti_list:
        m = fname.find( 'min' )
        low = int( fname[ m+3:m+7 ] )
        nx.append( low )

    bettis = []
    for b in betti_list:
        with open( fdir + b ) as fh:
            data = pkl.load( fh )
        bnums = data.values()[0][dim]
        bettis.append( bnums )
    bettis = np.array( bettis, dtype=np.int )

    if fig is None:
        fig = plt.figure(figsize=(8, 6)) 
    ax = fig.gca()
    ax.plot( bettis, color+shape+'-', ms=size, label=label )
    ax.set_xlabel( "Block number (height="+ height_ +")", fontsize=16 )
    ax.set_ylabel( r"$\beta_{"+str( dim )+"}$", fontsize=16 )
    ax.set_xticklabels( nx )
    ax.set_xlim( -1, len(bettis) )
    #ax.set_ylim( min(bettis)-1, max(bettis)+1 )

    ax.legend()

    return fig, (nx,bettis)

def make_combined_betti_figs( nx, data1, data2, height, dim ):
    """
    Comparison of capped and uncapped blocks. 

    Plot the homology of the blocks on top, the differences below.

    nx : Blocks base for x axis

    data1, data2 : 1D numpy arrays of betti numbers across blocks.

         data1 : uncapped

         data2 : capped (if these are switched, the difference plot
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

    print "data1", data1
    print "data2", data2

    print diffs

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
    # Build stacks
    #===========================
    if 0:

        betti_path = '/sciclone/data10/jberwald/CT_Firn_Samples/firn_blocks_June2013/'
        #betti_path = '/data/CT_Firn_Sample/firn_blocks/'

        #low = 3200
        #high = 3800

        h = 60
        # input from shell script
        low = int( sys.argv[1] )
        high = low + h

        print "height = ", h
        print low, high

        # do the homology computations
        B = firn_stacks( low, high, h, path, prefix, betti_path )

        # at each height, save the betti numbers of each stack
        savename = betti_path + prefix + 'min' + str(low) +\
            '_max'+str(high) + '_block' + str(h) + '.pkl'
        with open( savename, 'w' ) as fh:
            pkl.dump( B, fh )

        # print the homology so it shows up in output file
        print B
        print "Total time: ", time.time() - start


    #===========================
    # Build stacks with caps
    #===========================
    if 0:
        betti_path = '/sciclone/data10/jberwald/CT_Firn_Samples/firn_caps_June2013/'
        #betti_path = '/data/CT_Firn_Sample/firn_caps/'

        h = 40

        # input from shell script
        low = int( sys.argv[1] )
        high = low + h

        print "height = ", h
        print low, high

        B = firn_stacks( low, high, h, path, prefix, betti_path, cap=True )

        savename = betti_path + prefix + 'min' + str(low) +\
            '_max'+str(high) + '_block' + str(h) + '.pkl'
        with open( savename, 'w' ) as fh:
            pkl.dump( B, fh )

        print "Total time: ", time.time() - start

    if 1:

        betti_blocks = '/data/CT_Firn_Sample/firn_blocks/'
        betti_caps = '/data/CT_Firn_Sample/firn_caps/'
        dims = [0,1,2] #[0,1,2]
        stack_height = [40] #,60]

        for h in stack_height:
            print "stack height", h
            for d in dims:
                fig, b1 = make_betti_fig( h, d, betti_blocks, prefix, label='Uncapped' )
                fig, b2 = make_betti_fig( h, d, betti_caps, prefix, fig=fig, color='r', 
                                          shape='s', label='Capped' )

                fig.savefig( '/data/CT_Firn_Sample/figures/blocks_caps/'+\
                             'betti_h'+str(h)+'_d'+str(d)+'.pdf' )
                #fig.show()

        
                # now, take b1,b2 data and plot difference plots
                nx = b1[0]
                d1 = b1[1]
                d2 = b2[1]
                fig_comb = make_combined_betti_figs( nx, d1, d2, h, d )
                fig_comb.show()
