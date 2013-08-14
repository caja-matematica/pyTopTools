"""
pyImage.py

author: Jesse Berwald

opened: June 10, 2013

Library for manipulating image data with PIL and interfacing with
various libraries. Primarily aimed at homological analysis of images
at this point.
"""
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
import subprocess as sp
from itertools import izip

# local
import npy2perseus as n2p
import perseus_wrap as perseus


class PyImage( object ):
    """
    Class for dealing with BMP images.
    """
    def __init__( self, imagepath ):

        self.image = Image.open( imagepath )

    def bmp2array( self ):
        """
        """
        nx,ny = self.image.size # tuple (x,y)
        data = np.array( self.image.getdata(), dtype=np.int )
        data.resize( (ny,nx) ) # this must be reversed
        self.data = data
  
class ChompImage( PyImage ):

    def __init__( self, image ):
        """
        image : nd image array. typically, a 2d (single image) or 3d
        (sequence of images).
        """
        self.data = image

    def get_cubical_corners( self, val=255 ):
        """arr : binarized (thresholded) array.

        val : Calue to threshold on. Default==255. Based on greyscale
        images. In particular, firn images that are already
        thresholded to values of {0,255}. In the default case, 255 is
        'topological space', and 0 is not in the 'space' (i.e. it
        represents a hole!)

        Returns a list of indices suitable for writing a cubical complex.
        A cubical complex file consists of a list of many lines of the
        form ( n1, n2, ..., nd ) where d is the dimension of the complex
        and the ni's are the coordinates of the individual cube.

        Example: (0, 0, 1, 0) (10, 13, 2, 3) ... two cubes of a four
        dimensional cubical complex

        """
        w = np.where( self.data == val )
        self.corners = list( izip( *w ) )
    
    def cub2file( self, cubname ):
        """
        Convert an array of indices to chomp format, ( , , ). Write the resulting
        column of numbers to disk.
        """
        rows = map( lambda x: str(x)+'\n', map( tuple, iter( self.corners ) ) ) 
        with open( cubname, 'w' ) as fh:
            fh.writelines( rows )

            
    def run_chomp( self, fname ): #, savename ):
        """
        Call chomp to compute the betti numbers of the image in file fname.
        
        See http://chomp.rutgers.edu
        """
        try:      
            p = sp.check_output( ["chomp-rutgers", fname] )
        except:
            print "subprocess returned with command: ", cmd
        self.hom_output = p
        return p

    def extract_betti( self ):
        """
        chomp_out -- string output from chomp-rutgers
        """
        out = self.hom_output
        sout = out.split( '\n' )
        for line in sout:
            if line.startswith( 'Betti' ):
                # keep only the numbers
                # index=2 excludes the text
                betti_numbers = line.strip().split()[2:]
        self.betti = [ int( b ) for b in betti_numbers ]

class PerseusImage( PyImage ):
    
    def __init__( self, image ):
        """
        image : numpy nd array
        """
        self.data = image
        self.dim = image.ndim
    

    def write_perseus( self, persfile, scale=1, dtype='scubtop' ):
        """
        persfile : path to Perseus file format output
        """
        self.persfile = persfile
        if dtype == 'scubtop':
            n2p.write_scubtop_ND( self.data, persfile )
        elif dtype == 'cubtop':
            n2p.write_cubtop( self.data, persfile, scale )
            
    
    def run_perseus( self, persout, dtype='scubtop' ):
        """
        Run perseus on the data. Must run write_perseus() before calling.

        persout : Filename prefix for perseus output. 
        Eg., mypersfile --> mypersfile_d.txt,

        where d is the dimension. Note, perses add the '_*.txt'.
        """
        self.persout = persout
        perseus.perseus( self.persfile, self.persout, dtype=dtype )

    def draw_persdia( self, dim, scale=None, ms=1 ):
        """
        Draw the diagram for the dim'th persistent homology.

        ms : marker_scale
        """
        pdfile = self.persout + '_' + str(dim) + '.txt'
        if scale is None:
            fig = perseus.plot_diagram( pdfile )
        else:
            fig = perseus.plot_diagram_scaled( pdfile, fontsize=12, scale=scale,
                                               marker_scale=ms)
        return fig


#====================================

def cap_block( ims, cap_value=255, dtype=np.int, top=True, bottom=True ):
    """Cap a block of data with an array consisting of a single value.

    ims : list of image arrays or shape (nx,ny).

    cap_value : corresponds to the value of solid material in an
    image.

    dtype : default value for block array

    top : turn off the top cap

    bottom : turn off the bottom cap

    A capped block is returned |------------|, where | == array of
    constant value==cap_value

    Note: definitely duplicating ims in memory here!

    """
    if not (top==True or bottom==True):
        print "At one of 'top' or 'bottom' must be true in order to cap the block."
        return 

    # just in case we send in a list of 2D images
    ims = np.asarray( ims, dtype=dtype )

    # stack height == nx, (ny,nz)==image dims
    nx, ny, nz = ims.shape
    
    # update the dimensions
    if bottom:
        nx += 1
    if top:
        nx += 1

    block = np.empty( (nx, ny, nz), dtype=dtype )

    thecap = np.empty( (ny,nz), dtype=dtype )
    thecap.fill( cap_value )
        
    # place the cap on the block
    if bottom:
        block[0] = thecap
    if top:
        block[-1] = thecap

    # fill the rest of the block with arrays from ims. Duplicates all
    # of the block in memory!
    if top and bottom:
        block[1:-1] = ims 
    elif top:
        block[:-1] = ims
    else:
        block[1:] = ims

    return block
        

#################
# TESTING
#################

def test_PyImage():
    fname = 'sandbox/K09b-23-10-3230.bmp'
    B = PyImage( fname )
    B.bmp2array()
    return B

def test_ChompImage():
    fname = 'sandbox/K09b-23-10-3230.bmp'
    B = PyImage( fname )
    B.bmp2array()

    cb = ChompImage( B.data )
    cb.get_cubical_corners()
    cb.cub2file( 'sandbox/CUBTEST.cub' )
    cb.run_chomp( 'sandbox/CUBTEST.cub' )
    cb.extract_betti()
    return cb

def test_blocks():
    prefix = '/data/CT_Firn_Sample/output23-10-3/K09b-23-10-323'
    ims = []
    for i in range(1):
        im = PyImage( prefix + str(i) + '.bmp' )
        im.bmp2array()
        ims.append( im.data )

    ims = np.array(ims) 

    #block = build_block( ims )
    block = cap_block( ims )

    #print block.shape
    block = np.array( ims )
    return block, ims

def test_Chomp_blocks():
    bk = test_blocks()

    cb = ChompImage( bk )
    cb.get_cubical_corners()
    cb.cub2file( 'sandbox/CUBTEST.cub' )
    cb.run_chomp( 'sandbox/CUBTEST.cub' )
    cb.extract_betti()
    return cb

def test_block_small( top=True, bottom=True, test_num=0 ):
    ims = [ np.zeros( (3,3), dtype=np.int )
            for i in range( 5 ) ]

    for im in ims: 
        im.fill( 255 )
    
    # Now create a hole in the center (val!=255 will do it)
    if test_num == 0:
        for im in ims:
            im[1,1] = 0
    # L-shaped divot removed from 2D images
    elif test_num == 1:
        for im in ims:
            im[1,1] = im[1,2] = im[2,1] = 0

    ims = np.array( ims )

    if top or bottom:
        c = cap_block( ims, top=top, bottom=bottom )
    else:
        c = ims        
    cb = ChompImage( c )
    cb.get_cubical_corners()
    cb.cub2file( 'sandbox/SMALLTEST.cub' )
    cb.run_chomp( 'sandbox/SMALLTEST.cub' )
    cb.extract_betti()
    return cb, c
    
def test_Chomp_image_cap():
    """Stack four images and cap them
    """
    prefix = '/data/CT_Firn_Sample/output23-10-3/K09b-23-10-323'
    ims = []
    for i in range(1):
        im = PyImage( prefix + str(i) + '.bmp' )
        im.bmp2array()
        ims.append( im.data ) 

    #cap = cap_block( ims )
    cap = np.array( ims )

    cb = ChompImage( cap )
    cb.get_cubical_corners()
    cb.cub2file( 'sandbox/CUBTEST.cub' )
    cb.run_chomp( 'sandbox/CUBTEST.cub' )
    cb.extract_betti()
    return cb

def test_cap():
    """
    Test that a 3x3x5 block, with a 'hole' in the center 3 arrays
    gives \beta_0=1, \beta_1=0, \beta_2=1. [passed]
    """
    ims = [ np.zeros( (3,3), dtype=np.int )
            for i in range( 3 ) ]
    for im in ims: 
        im.fill( 255 )    
    for im in ims:
        im[1,1] = 1
    c = cap_block( ims, cap_value=255 )
    cb = ChompImage( c )
    cb.get_cubical_corners()
    cb.cub2file( 'sandbox/CAPTEST.cub' )
    cb.run_chomp( 'sandbox/CAPTEST.cub' )
    cb.extract_betti()
    return cb, c

def test_persistence():
    """
    """
    from pyTopTools.rbc_analysis import gaussian as G
    matfile = './rbc_analysis/data_for_figs/gauss_smooth.npy'
    data = np.load( matfile )

    # set values below 1 to 0 (in-place)
    G.clip_below( data, 1 )

    persfile = './sandbox/PERSTEST.txt'
    persout = './sandbox/PERSOUT'
    
    P = PerseusImage( data )
    P.write_perseus( persfile )
    P.run_perseus( persout )
    return P
    

def test_pers_stack():
    """
    Stack a number of Gaussians
    """
    from pyTopTools.rbc_analysis import gaussian as G
    matfile = './rbc_analysis/data_for_figs/gauss_smooth.npy'
    data = np.load( matfile )

    # set values below 1 to 0 (in-place)
    # for scubtop
    G.clip_below( data, 1 )

    nx,ny = data.shape
    block = np.zeros( (nx,ny,5) )
    for i in range(5):
        block[:,:,i] = data

    persfile = './sandbox/PERSTEST_3D.txt'
    persout = './sandbox/PERSOUT_3D'

    P = PerseusImage( block )
    P.write_perseus( persfile )
    P.run_perseus( persout )
    return P

def test_pers_small( cap=False ):
    """A block with two tubes. So

    H_0 = {1..-1}
    H_1 = {1..10}
    H_2 = 0

    If cap==True, then

    H_0 = {1..-1}
    H_1 = 0
    H_2 = {1..-1}

    Note, cap_val=1, since Perseus interprets vals <= 0 as never
    existing in the space.

    """
    a = np.ones( (10,10) )
    a[3,3] = a[7,7] = 10
    
    b = [ a for i in range(3) ]
    
    if cap:
        c = cap_block( b, cap_value=1 )
    else:
        c = b
    
    persfile = './sandbox/PERSTEST.txt'
    persout = './sandbox/PERSOUT'
    
    P = PerseusImage( c )
    P.write_perseus( persfile )
    P.run_perseus( persout )
    return P
    
    

if __name__ == "__main__":

    # p = test_PyImage()
    #c = test_ChompImage()
    #b = test_blocks()
    #cc = test_Chomp_blocks()
    #ccap = test_Chomp_image_cap()
    

    test = 1
    # cap only top, no bottom
    #cct, ccimsB = test_block_small( bottom=False, test_num=test )

    # cap bottom, no top
    #ccb, ccimsT = test_block_small( top=False, test_num=test )

    # cap both ends
    #cc, ccims = test_block_small( test_num=test )

    # no cap
    #cno, cnoims = test_block_small( top=False, bottom=False, test_num=test )

    #cp, ims = test_cap()

    pers = test_persistence()
    #pers3 = test_pers_stack()
    #psmall = test_pers_small( cap=True )
