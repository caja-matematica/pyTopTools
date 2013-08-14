"""
Compute Wasserstein distance for various lags.
"""
import sys
import numpy as np
import subprocess as sp
import cPickle as pkl
from pylab import norm

def write_robust_gens( arr1, arr2, out1, out2, fmt='%g' ):
    """Write robust generator array to file, eg., peristence
    diagram. Used in conjunction with compute_wasserstein() to write
    two lagged persistence diagram to file for computing the distance
    between them.

    arr* : the persistence diagram arrays. 

    out* : output arr* to out*. out* must contain the frame number to
    guarantee unique file names.

    """
    np.savetxt( out1, arr1, fmt=fmt )
    np.savetxt( out2, arr2, fmt=fmt )
    
def compute_wasserstein( pd1, pd2, p=2 ):
    """
    pd1, pd2 : persistence diagrams to compute distance between.

    p : Compute the p-Wasserstein distance
    """
    out = sp.check_output( ["wasserstein", pd1, pd2] )
    return out

def lag_vector( vector_path, p=2 ):
    """
    vector_path : path to lag vector on disk, 

    eg., /sciclone/data10/jberwald/RBC/cells/persout/old_8_pdist_lag1.npy

    p : int or 'inf'
    """
    vec = np.load( vector_path )
    if p == 'inf':
        vecnorm = norm( vec, ord=np.inf )
    else:
        vecnorm = norm( vec, ord=p )
    return vecnorm
        
def frame_max( frame, bndy ):
    """Returns vector of max of frame from a time series of cell, with the
    boundary removed.

    frame : single frame from cell time series

    bndy : cell boundary matrix

    """
    cell = frame.ravel()
    cell_clipped = bndy * cell
    return cell_clipped.max() 

def normalize_gens( pd, frame_max ):
    """pd : persistence diagram, an n x 2 numpy array
    
    frame_max : maximum value for the frame giving associated to the
    pd.

    """
    pd = np.asarray( pd, dtype=np.float )
    return pd / frame_max

def trim_longest( pd ):
    """Remove the row of the pd containing the 'infinite' generator.
    """
    # lifespans
    diff = np.diff( pd )
    max_life = diff.argmax()
    # delete row=max_life 
    new_pd = np.delete( pd, (max_life), axis=0 )
    return new_pd
    

def main( cell, prefix, p=2, lag=1, out_prefix='/sciclone/scr01/jberwald/' ):
    """Routine to compute the Wasserstein distance between cell frames
    with a lag of 'lag'.
    """
    if cell.startswith('old'):
        suffix = 'old_robust.pkl'
    else: 
        suffix = 'new_robust.pkl'
    with open( prefix + suffix ) as fh:
        all_robust = pkl.load( fh )

    # we're analyzing only frames in 'cell'
    diagrams = all_robust[ cell ]

    # save prefix for txt output
    out_prefix += cell
    out1 = out_prefix + '_1_'+str( lag )
    out2 = out_prefix + '_2_'+str( lag )

    frames = diagrams.keys()
    distances = []
    for frame in frames[:-lag]:
        pd1 = diagrams[ frame ]
        pd2 = diagrams[ frame+lag ]
        
        # delete the longest (inf) generator
        pd1 = trim_longest( pd1 )
        pd2 = trim_longest( pd2 )

        # just store each pair and overwrite at next step
        write_robust_gens( pd1, pd2, out1, out2 )
        pdist = compute_wasserstein( out1, out2 )
        distances.append( int( pdist.strip() ) )

    D = np.asarray( distances )
    np.save( prefix + cell + '_pdist_lag'+str(lag), D )
    return D


def main_normed( cell, prefix, p=2, lag=1, 
                 out_prefix='/sciclone/scr00/jberwald/',
                 cell_suffix='-concatenated-ASCII',
                 remove_longest=False ):
    """Normalize the generators in each frame by the maximum value of each
    frame. Then compute the Wasserstein distance between frames. Uses
    only robust generators with lifespans > 40.

    cell : single cell name (to be appended to out_prefix)

    prefix : path to directory containing robust generator information

    p : p-norm

    lag : distance between compared frames. Eg., 
    
    d( frame_i, frame_{i+lag} )

    """
    if cell.startswith('old'):
        suffix = 'old_robust.pkl'
    else: 
        suffix = 'new_robust.pkl'
    with open( prefix + suffix ) as fh:
        all_robust = pkl.load( fh )

    # save prefix for txt output
    out_prefix += cell
    out1 = out_prefix + '_1_'+str( lag )
    out2 = out_prefix + '_2_'+str( lag )

    # we're analyzing only frames in 'cell'
    diagrams = all_robust[ cell ]

    frames = diagrams.keys()
    distances = []
    for frame in frames[:-lag]:
        # these include the longest (inf) generators
        pd1 = diagrams[ frame ]
        pd2 = diagrams[ frame+lag ]

        # normalize by the max pixel value (which should be the death
        # of the inf generator). grab this before we remove it below
        pd1max = pd1.max()
        pd2max = pd2.max()

        if remove_longest:
            # returns diagram with longest lifespan generator removed
            pd1 = trim_longest( pd1 )
            pd2 = trim_longest( pd2 )

        # now do the normalization
        pd1 = normalize_gens( pd1, pd1max )
        pd2 = normalize_gens( pd2, pd2max )

        # just store each pair and overwrite at next step
        write_robust_gens( pd1, pd2, out1, out2 )
        pdist = compute_wasserstein( out1, out2 )
        distances.append( float( pdist.strip() ) )

    D = np.asarray( distances )
    np.save( prefix + cell + '_pdist_lag'+str(lag)+'_norm', D )
    return D

if __name__ == "__main__":

    # usage: python compute_wasserstein.py cell_name lag
    cell_name = sys.argv[1]
    lag = int( sys.argv[2] )

    prefix = '/sciclone/data10/jberwald/RBC/cells/persout/'
    dists = main_normed( cell_name, prefix, lag=lag, remove_longest=True )
    dists = main( cell_name, prefix, lag=lag )
