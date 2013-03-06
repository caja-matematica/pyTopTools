import numpy as np
import matplotlib.pyplot as plt
from pyRBC import rbc_postprocess as rpost
from pyRBC import rbc_npy2Perseus as rn
from pyRBC import rbc_perseus as rp
from pyRBC import rbc_perseus as rpers


class Frame( object ):
    """
    """
    def __init__( self, the_frame ):
        """
        the_frame : a single, resized (n x m) frame (numpy ndarray).
        """
        self.frame = the_frame

    def write_frame_npy( self, fnum ):
        """
        Write frame to disk in .npy format.
        """
        np.save( self.fname + under + str( frame_num ) + '.npy', self.frame )
    
    def write_perseus( self, frame, outname ):
        rn.write_sparse_array( frame, outname )

    def run_perseus( self, pers_in, pers_out ):
        rpers.perseus( pers_in, pers_out )

    def draw_persdia( self, fname, **args ):
        rpost.plot_diagram_std( fname, **args )

    def imshow( self ):
        fig = plt.figure()
        ax = fig.gca()
        ax.imshow( self.frame )


class Cell( Frame ):

    def __init__( self, fname, bndy, skiprows=0, shape=(203,198), verbose=False ):
        """
        Load a complete, concatenated cell text file.

        fname : path to file containing concatenated cell frames. 

        skiprows : skip the first 'skiprows'.

        shape : shape for resizing each frame. default=(203,198)

        square : dimension for reshaping boundary to and cell frames to a
        square matrix. necessary for lining them up, since even though
        boundary is centered correctly over cell, this only occurs when
        shapes of matrices _do not_ line up.
        """
        self.fname = fname
        self.shape = shape
        self.verbose = verbose
        
        self.all_frames = np.loadtxt( fname, skiprows=skiprows, dtype=np.uint )
        if self.verbose:
            print "Loaded frames"
            
        self.frames = []
        for frame in self.all_frames:
            frame.resize( shape )
            self.frames.append( Frame( frame ) )

        self.boundary = np.loadtxt( bndy )#, dtype=np.uint )
        self.boundary.resize( (shape[1],shape[0]) )
        
        if self.verbose:
            print "Done initializing frames and boundary"


if __name__ == "__main__":
    
    # TESTING
    
    prefix = '/data/jberwald/rbc/concat_files/old_6/'
    cell_name = 'old_6-concatenated-ASCII'
    bnd = 'boundary_April_old6-col'
    persname = prefix + 'noise/old_6_pers.txt'
    persout = prefix + 'noise/old_6_pers'

    C = Cell( prefix + cell_name, prefix + bnd, skiprows=4990 )
