import numpy as np
from pyRBC import rbc_npy2Perseus as rn
from pyRBC import rbc_perseus as rp
from pyRBC import rbc_postprocess as rpost
from pyRBC import rbc_perseus as rpers


"""
Module specifically written for cell old_6. Reshaping of boundary and
cell specifically works for this cell. (Untested on others.)
"""

under = '_'

class Cell( object ):

    def __init__( self,  fname, bndy, skiprows=0, shape=(203,198), verbose=False ):
        """
        Load a complete, concatenated cell text file.

        fname : path to file. 

        skiprows : skip the first 'skiprows'.

        shape : shape for resizing each frame. default=(203,198)

        square : dimension for reshaping boundary to and cell frames to a
        square matrix. necessary for lining them up, since even though
        boundary is centered correctly over cell, this only occurs when
        shapes of matrices _do not_ line up.
        """
        self.fname = fname
        self.shape = shape
        if shape[0] > shape[1]:
            self.square = shape[1]
        else:
            self.SQUARE = shape[0]
        self.verbose = verbose
        
        self.all_frames = np.loadtxt( fname, skiprows=skiprows, dtype=np.uint )
        if self.verbose:
            print "Loaded frames"
            
        self.frames = []
        for frame in self.all_frames:
            frame.resize( shape )
            self.frames.append( frame[:self.square,:self.square] )

        B = np.loadtxt( bndy )#, dtype=np.uint )
        B.resize( (shape[1],shape[0]) )
        self.boundary = B[:self.square,:self.square]

        if self.verbose:
            print "Done initializing frames and boundary"

    def _reshape( self, A ):
        return A.reshape( (self.square,self.square) )


    def crop_cell( self, fnum, threshold=800 ):
        """
        Set all elements within cell boundary to zero.
        """
        # first clip the portion of the cell inside the halo by
        # creating an inverted boundary matrix.
        invert_bnd = np.zeros_like( self.boundary, dtype=np.uint )
        out_halo = np.where( self.boundary == 0 )
        invert_bnd[ out_halo ] = 1

        # now clip the cell image. this leaves the halo...
        invC = self.frames[ fnum ] * invert_bnd

        # trim the halo uing 'threshold'
        thresh = np.where( invC > threshold )
        invC[ thresh ] = 0
        return np.asarray( invC, dtype=np.uint )
        
    def write_frame_npy( self, fnum ):
        """
        frames : list of frames (numpy matrices) to save to disk
        """
        np.save( self.fname + under + str( frame_num ) + '.npy', self.frames[ fnum ] )
    
    def write_perseus( self, frame, outname ):
        if type( frame ) is int:
            rn.write_sparse_array( self.frames[ fnum ], outname )
        else:
            rn.write_sparse_array( frame, outname )

    def run_perseus( self, pers_in, pers_out ):
        rpers.perseus( pers_in, pers_out )

    def draw_persdia( self, fname, **args ):
        rpost.plot_diagram_std( fname, **args )


if __name__ == "__main__":

    prefix = '/data/jberwald/rbc/concat_files/old_6/'
    cell_name = 'old_6-concatenated-ASCII'
    bnd = 'boundary_April_old6-col'
    persname = prefix + 'noise/old_6_pers.txt'
    persout = prefix + 'noise/old_6_pers'

    C = Cell( prefix + cell_name, prefix + bnd, skiprows=4990 )
    clip = C.crop_cell( 0 ) # crop the first frame
    C.write_perseus( clip, persname )
    C.run_perseus( persname, persout )
    
    
