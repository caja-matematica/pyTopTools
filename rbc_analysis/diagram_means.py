import numpy as np
import cPickle as pkl
import matplotlib.pyplot as plt

# local
from pyTopTools import persistence_tools as pt
from pyTopTools import persistence_diagrams as pd

def sample_frames( robust, prefix, trial, sample_size=100 ):
    """
    robust : dictionary collection of all diagrams containing robust
    generators (lifespan > 40)

    prefix : directory in which to save samples (in matlab format for
    diagram mean computations)

    trial : sample number

    Optional:
    --------

    sample_size : number of sample frames to choose.
    """
    
    # subsample frames
    for cell in robust.keys():
        frames = robust[cell].keys()

        # random sample of frame to find the mean for
        xmin = min( frames ) # should be 0
        xmax = max( frames ) # a few cell have < 5000 frames
        samples = np.random.random_integers( xmin, xmax, sample_size )

        # create list of diagrams to pass to cellarray converter
        dia_list = [ robust[cell][ i ] for i in samples ]

        # create name of .mat file
        sample_output = cell + '_dia_samples_' + str(trial)
        outname = prefix + sample_output
        pt.diagrams2cellarray( dia_list, outname, chop_inf=True, mat_type=np.float )

def mean_distances( mean_sample, dim=1, matname='means' ):
    """
    cluster cells based on mean diagrams
    """
    if not mean_sample.endswith('.mat'):
        mean_sample += '.mat'
    means = pt.cellarray2diagrams( mean_sample, matname=matname )
    
    distances = np.zeros( ( len(means), len(means) ) )
    cells = means.keys()

    cell_hash = {}
    computed = []
    diagrams = []
   
    for i,cell in enumerate( cells ):
        # keep track of which cell is in which row
        cell_hash[i] = cell
        D = pd.Diagram( means[ cell ], dim=dim )
        diagrams.append( D )
        for j,other in enumerate( cells ):
            if i==j:
                continue
            if j > i and ([i,j] in computed):
                distances[j,i] = distances[i,j]
                computed.append( [j,i] )
                continue

            D2 = pd.Diagram( means[ other ] )
            d =  D.compute_wasserstein_distance( D2 )
            # d =  D.compute_wasserstein_distance( D2 )

            distances[i,j] = d # D.compute_wasserstein_distance( D2 )
            computed.append( [i,j] )
    return distances, cell_hash, diagrams

def convert_to_compressed_distance_array( D ):
    """
    Convert square distance array D into a m(m-1)/2 vector for Mapper.

    See mapper doc: http://math.stanford.edu/~muellner/mapper/input.html

    """
    # get indices of upper triangular matrix
    n,n = D.shape # assume square!
    idx = np.triu_indices( n, k=1 ) # diagonal should be zero
    return D[idx]


def test_distance( distmat, dimcap=3, dtype='distmat', **kwargs ):
    """
    """
    from scipy.spatial.distance import squareform
    from pyPersistence import filter_data as FD

    persfile = './data/DISTTEST.txt'
    persout = './data/DISTOUT'

    # correlation matrix
    data = np.loadtxt( distmat )
    data = squareform( data )

    # initialize PP object
    P = FD.PerseusProcessor( data )
    
    print "Testing write_perseus() for dtype=", dtype
    print "-----"
    P.write_perseus( persfile, dtype=dtype, dimcap=dimcap,
                     **kwargs )
    
    print "Calling Perseus on", persfile
    P.run_perseus( persout, dtype=dtype )
    print "Success!"    
    print "-----"

    return P


def plot_diagram_mean_sample( dias, cell_hash ):
    """
    dias : list of Diagram objects

    cell_hash : list index --> cell name mapping
    """
    fig = plt.figure()

    for ii, cellname in cell_hash.items():
        D = dias[ii]
        if 'old' in cellname:
            c = 'r'
        else:
            c = 'b'
        D.draw_diagram( fig=fig, color=c, marker_scale=8 )
    
    fig.show()
    return fig

if __name__ == "__main__":

    sample_size = 100
    trials = np.arange( 500 )

    #prefix = '/Users/jberwald/github/local/caja-matematica/pyTopTools/rbc_analysis/data/'
    prefix = '/ima/imausr1/jberwald/data/rbc/mean_dias/'
    means_prefix = prefix 
    #'/ima/imausr1/jberwald/Dropbox/Work/Topology/wassMean/means/'

    # with open( prefix + 'robust_eps40.pkl' ) as fh:
    #     robust = pkl.load( fh )

    for t in trials:
        print "computing for trial", t

        #sample_frames( robust, prefix, t )

        fname = 'mean_sample_'
        D, cell_hash, diagrams = mean_distances( means_prefix +\
                                                     fname + str(t)+'.mat',
                                                 matname='cellMeans' )
        Dtri = convert_to_compressed_distance_array( D )
        np.savetxt( means_prefix+'mean_dia_distances_' + str(t) +'.txt', Dtri )
        with open( means_prefix+'mean_dias_'+str(t)+'.pkl', 'w' ) as fh:
            pkl.dump( diagrams, fh )
        with open( means_prefix+'cellhash_dias_'+str(t)+'.pkl', 'w' ) as fh:
            pkl.dump( cell_hash, fh )
   
    #D, cell_hash, diagrams = mean_distances( means_prefix + 'mean_sample_1.mat' )
    
    
