import numpy as np
import cPickle as pkl
import tempfile
import time
import matplotlib.pyplot as plt

from pyTopTools import persistence_diagrams as pd
from pyTopTools import tools

prefix = '/ima/imausr1/jberwald/data/rbc/'


def compute_mutual_distances():
    print "Loading files..."

    with open( prefix + 'old_robust.pkl' ) as fh:
        oldcells = pkl.load( fh )
    with open( prefix + 'new_robust.pkl' ) as fh:
        newcells = pkl.load( fh )
    print "...done.\n"

    newnames = newcells.keys()
    oldnames = oldcells.keys()
    allnames = newnames + oldnames

    mutual_dists = dict()

    #allnames = ['old_8', 'new_50125']

    for cell1 in allnames:
        try:
            dia_list1 = newcells[ cell1 ]
        except KeyError:
            dia_list1 = oldcells[ cell1 ]
        for cell2 in allnames:
            if cell1 == cell2: continue
            try:
                dia_list2 = newcells[ cell2 ]
            except KeyError:
                dia_list2 = oldcells[ cell2 ] 

            print "Computing diagram distance along time series for cells:"
            print cell1
            print cell2 
            print ""

            # align time since some time series are a little short
            N = min( len(dia_list1), len( dia_list2 ) )
            time_arr = np.arange( 0, N, 10 )
            # the 'inner' loop
            W = []
            for i in time_arr:
                D1 = pd.Diagram( dia_list1[i], dim=1 )
                D2 = pd.Diagram( dia_list2[i], dim=1 )
                dist = D1.compute_wasserstein_distance( D2 )
                W.append( int(dist) )

            mutual_dists[ ( cell1, cell2 ) ] = W    
    return mutual_dists

def mutual_sums( fname ):
    """
    fname : path to pkl file containing mutual distance between time
    series of persistence diagrams.
    """
    with open( fname ) as fh:
        M = pkl.load( fh )

    # array to hold variables
    dt = np.dtype([('cell1', np.str_, 32), 
                   ('cell2', np.str_, 32),
                   ('dsum', np.float64, (1,) )])
    D = np.empty( (26**2, 3), dtype=dt )
    # keys == ('cell1', cell2')
    for cells, dia_dists in M.items:
        
        total_distance = sum( dia_dists )
        

def idx_of_inf_gen( dia ):
    """
    n x 2 diagram array
    """
    d = np.diff( dia )
    return d.argmax() 

def compute_lifespan_sums( noise_factor=40 ):
    """
    """
    print "Loading files..."

    with open( prefix + 'old_robust.pkl' ) as fh:
        oldcells = pkl.load( fh )
    with open( prefix + 'new_robust.pkl' ) as fh:
        newcells = pkl.load( fh )
    print "...done.\n"

    newnames = newcells.keys()
    oldnames = oldcells.keys()
    allnames = newnames + oldnames

    total_persistence = {}

    for cell in allnames:
        try:
            frames = oldcells[cell]
        except KeyError:
            frames = newcells[cell]
        cell_sums = []
        for dia in frames.values():
            # remove the infinite generator 
            idx = idx_of_inf_gen( dia )
            robust = np.delete( dia, idx, axis=0 )
            robust_diff = np.diff( robust ) - noise_factor
            cell_sums.append( robust_diff.sum() )
        total_persistence[ cell ] = cell_sums
    return total_persistence

def compute_mse( data, savename ):
    """
    """
    fh = tempfile.NamedTemporaryFile( delete=False, prefix='/scratch/tmp' )
    np.savetxt( fh, data )
    fh.close()
    tools.run_mse( fh.name )
    time.sleep(1) # make sure the file has been written before we try to read it
    d = tools.mse_converter( fh.name + '.mse' )
    np.savetxt( savename, d )
    return d 

def plot_mse( ):
    pass

if __name__ == "__main__":
    
    # M = compute_mutual_distances()
    # with open( 'mutual_diagram_distances.pkl', 'w' ) as fh:
    #     pkl.dump( M, fh )

    print "Computing lifespan sums..."
    T = compute_lifespan_sums()
    LS = {}
    # for k,v in T.items():
    #     print "Computing MSE for lifespans of ", k
    #     print ""
    #     s = compute_mse( v, './data/'+k+'.mse' )
    #     LS[ k ] = s 

    # fig = plt.figure()
    # ax = fig.gca()
    # for k,v in LS.items():
    #     v = v.T
    #     if k.startswith( 'old' ):
    #         ax.plot( v[0], v[1], 'r-o', ms=6, lw=2 )
    #     else:
    #         ax.plot( v[0], v[1], 'b-o', ms=6, lw=2 )
    # ax.set_xlabel( 'Scale Factor', fontsize=16 )
    # ax.set_ylabel( 'Sample Entropy', fontsize=16 )
    # ax.set_title( 'MSE for lifespan sums over all cells (old=red, new=blue)',
    #               fontsize=16 )
    # plt.show()
