import numpy as np
import cPickle as pkl

# local
from pyTopTools import persistence_tools as pt
from pyTopTools import persistence_diagrams as pd

sample_size = 100
trial = 1

prefix = '/Users/jberwald/github/local/caja-matematica/pytoptools/rbc_analysis/data/'

with open( prefix + 'robust_eps40.pkl' ) as fh:
    robust = pkl.load( fh )

# subsample frames
for cell in robust.keys():
    frames = robust[cell].keys()
    # random sample of frame to find the mean for
    samples = np.random.choice( frames, sample_size )
    dia_list = [ robust[cell][ i ] for i in samples ]

    # create name of .mat file
    sample_output = cell + '_dia_samples_' + str(trial)
    outname = prefix + sample_output
    pt.diagrams2cellarray( dia_list, outname, chop_inf=True, mat_type=np.float )
