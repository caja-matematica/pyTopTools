from pyTopTools import chomp_image as C
from pyTopTools import pyImage as P
import numpy as np
import matplotlib.pyplot as plt
import cPickle as pkl
import pp, time


#chomp_path = '/sciclone/data10/jberwald/CT_Firn_Samples/chomp_files2/'
#chomp_path = './_cubfile.cub'

def firn_stack_betti_numbers( ims, cubfile='/var/tmp/_cubfile.cub' ):
    """
    Create a single stack of images and compute the homology.
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

def build_firn_stacks( bottom, top, height, step=None, cap=False ):
    """
    Build a sequence of image stacks:

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
    for base in range( bottom,  top, step ):
        # create a sequence of images to stack
        ims = []
        # grab each image from base --> base+height
        for x in range( height ):
            image_num = base + x
            im = P.PyImage( path + prefix + str( image_num ) + '.bmp' )
            im.bmp2array()
            ims.append( im.data )

        if cap:
            block = P.cap_block( ims )
        else:
            block = P.build_block( ims )
        # returns a ChompImage object
        b = firn_stack_betti_numbers( block )
        key = ( base, base+x )
        betti[ key ] = b.betti
        
    return betti


###################################
#
# Run experiments below
#
###################################

start = time.time()

#path = '/data/CT_Firn_Sample/output23-10-3/'
#chomp_path = '/data/CT_Firn_Sample/debug/'
#prefix = 'K09b-23-10-'

path = '/sciclone/data10/jberwald/CT_Firn_Samples/output23-10-3/'
prefix = 'K09b-23-10-'
    
if 1:

    betti_path = '/sciclone/data10/jberwald/CT_Firn_Samples/firn_blocks_June2013/'
    
    low = 3200
    high = 3800
    height = [ 40, 60 ]

    for h in height:
        for level in range( low, high, h ):
            print "Building stack of images... "
            print "  level ", str( level )
            print "  height ", str( h )
            B = build_firn_stacks( low, high, h )
        savename = betti_path + prefix + 'base' + str(low) +\
            '_top'+str(high) + '_h' + str(h) + '.pkl'
        with open( savename, 'w' ) as fh:
            pkl.dump( B, fh )

    print "Total time: ", time.time() - start

if 0:

    betti_path = '/sciclone/data10/jberwald/CT_Firn_Samples/firn_caps_June2013/'
    
    low = 3200
    high = 3800
    height = [ 40, 60 ]

    for h in height:
        for level in range( low, high, h ):
            print "Building stack of images... "
            print "  level ", str( level )
            print "  height ", str( h )
            B = build_firn_stacks( low, high, h, cap=True )
        savename = betti_path + prefix + 'base' + str(low) +\
            '_top'+str(high) + '_h' + str(h) + '.pkl'
        with open( savename, 'w' ) as fh:
            pkl.dump( B, fh )

    print "Total time: ", time.time() - start

if 0:

    #for height in stack_height:
    height = 10
    dim = 1
    bettis = []
    for base in range( 3310,  3330, height ):
        betti_file = chomp_path + prefix[:-1] + \
            '_b' + str( base ) + \
            '_h' + str( height ) + \
            '.betti'
        print betti_file
        bnums = np.loadtxt( betti_file, dtype=np.uint8 )
        print bnums
        bettis.append( bnums[dim][1] )

    fig = plt.figure()
    ax = fig.gca()
    ax.plot( bettis, 'bo' )
    ax.set_xlabel( "Block number (height="+str(height)+")" )
    ax.set_ylabel( r"$\beta_{"+str( dim )+"}$" )
    ax.set_xlim( -1, len(bettis) )
    ax.set_ylim( min(bettis)-1, max(bettis)+1 )
    plt.show()
