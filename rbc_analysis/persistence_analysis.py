"""
Run Perseus on each from of a single cell file. 

author : Jesse Berwald

opened : June 28, 2013
"""
import numpy as np
import cPickle as pkl
import os, sys
from pyTopTools.rbc_analysis import rbc_perseus as rp
from pyTopTools import npy2perseus as n2p


def run_perseus( persin, persout, dtype='scubtop' ):
    """Run Perseus sparse cubical toplex analysis on a single frame.
    """
    rp.perseus( persin, persout, dtype )

def frame2scubtop( arr, persin ):
    """
    arr : single frame from concatenated cell file

    persin : name of persistence file to read by Perseus
    """
    n2p.write_scubtop_ND( arr, persin )

def read_persistence_diagram( fname ):
    return np.loadtxt( fname, dtype=np.int )

def compute_persistence( frame, persin, persout ):
    """
    """
    frame2scubtop( frame, persin )
    run_perseus( persin, persout )


if __name__ == "__main__":

    # for debugging
    skiprows = 0

    root_dir = '/sciclone/data10/jberwald/RBC/cells/'

    # prepare path to cell. cell name/number
    cell_prefix = sys.argv[1]
    cell_suffix = '-concatenated-ASCII'

    if cell_prefix.startswith( 'new' ):
        cellpath = root_dir + 'new/' + cell_prefix + cell_suffix
    else:
        cellpath = root_dir + 'old/' + cell_prefix + cell_suffix

    # prepare path to boundary file
    bndpath = root_dir + 'boundaries/boundary_' + cell_prefix
    boundary = np.loadtxt( bndpath )

    # this tell use how to reshape each frame from the cell file
    shape = boundary.shape

    # load the entire cell block into memory
    frames = np.loadtxt( cellpath, skiprows=skiprows )
    frames = [ x.reshape( shape ) for x in frames ]

    # trim the boundary off of the cell image. these are the data we will
    # work with
    cell = [ boundary * x for x in frames ]

    pers_root = root_dir + 'persout/' + cell_prefix
    PD = {}
    for i,frame in enumerate( cell ):
        pin = pers_root + '/'+ cell_prefix + '_frame'+str(i) + '_pers.txt'
        pout = pers_root + '/' + cell_prefix + '_frame'+str(i) + '_pers' # perseus will append '_*.txt'
        #compute_persistence( frame, pin, pout )
        PD[i] = read_persistence_diagram( pout + '_1.txt' )
        
    with open( pers_root + '_gens.pkl', 'w' ) as fh:
        pkl.dump( PD, fh )
