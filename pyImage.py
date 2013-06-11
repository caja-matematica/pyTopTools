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
        data = np.array( self.image.getdata(), dtype=np.uint8 )
        data.resize( (ny,nx) ) # this must be reversed
        self.data = data

  
class ChompImage( PyImage ):

    def __init__( self, image ):
        """
        image : nd image array. typically, a 2d (single image) or 3d
        (sequence of images).
        """
        self.data = image

    def get_cubical_corners( self, val=0 ):
        """
        arr : binarized (thresholded) array.

        Returns a list of indices suitable for writing a cubical complex.
        A cubical complex file consists of a list of many lines of the
        form ( n1, n2, ..., nd ) where d is the dimension of the complex
        and the ni's are the coordinates of the individual cube.

        Example: (0, 0, 1, 0) (10, 13, 2, 3) ... two cubes of a four
        dimensional cubical complex

        Note: This only works for 2D complexes
        """
        w = np.where( self.data == 0 )
        self.corners = zip( *w )
    
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



def cap_block( ims, cap_value=0 ):
    """
    Cap a block of data with an array consisting of a single value.

    ims : list of image arrays or shape (nx,ny).

    A capped block is returned |------------|, where | == array of
    constant value==cap_value
    """
    nx, ny = ims[0].shape
    block = np.zeros( ( nx,ny,len(ims)+2 ), dtype=np.uint8 )
    if cap_value == 0:
        thecap = np.zeros( (nx,ny), dtype=np.uint8 )
    else:
        thecap = np.emtpy( (nx,ny), dtype=np.uint8 )
    block[:,:,0] = block[:,:,-1] = thecap
    for i,im in enumerate( ims ):
        block[:,:,i+1] = im
    return block
                   

def build_block( ims ):
    """
    Frames : list of PyImage objects

    return n+1 dimensional block, where n=dim of each frame
    """
    nx, ny = ims[0].shape
    block = np.zeros( ( nx,ny,len(ims) ), dtype=np.uint8 )
    for i,im in enumerate( ims ):
        block[:,:,i] = im
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
    C = ChompImage( B.data )
    return C

def test_blocks():
    prefix = '/data/CT_Firn_Sample/output23-10-3/K09b-23-10-323'
    ims = []
    for i in range(4):
        im = PyImage( prefix + str(i) + '.bmp' )
        im.bmp2array()
        ims.append( im.data )
    block = build_block( ims )
    return block

def test_Chomp_blocks():
    bk = test_blocks()
    cb = ChompImage( bk )
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
    ims = [ np.zeros( (3,3), dtype=np.uint8 )
            for i in range( 3 ) ]
    for im in ims:
        im[1,1] = 1
    c = cap_block( ims )
    cb = ChompImage( c )
    cb.get_cubical_corners()
    cb.cub2file( 'sandbox/CAPTEST.cub' )
    cb.run_chomp( 'sandbox/CAPTEST.cub' )
    cb.extract_betti()
    return cb
    

if __name__ == "__main__":

    # p = test_PyImage()
    # c = test_ChompImage()
    # #b = test_blocks()
    # cc = test_Chomp_blocks()
    cp = test_cap()

