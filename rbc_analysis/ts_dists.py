import numpy as np
import cPickle as pkl
from pyTopTools import persistence_diagrams as pd

prefix = '/ima/imausr1/jberwald/data/rbc/'

with open( prefix + 'old_robust.pkl' ) as fh:
    oldcells = pkl.load( fh )
with open( prefix + 'new_robust.pkl' ) as fh:
    newcells = pkl.load( fh )

newnames = newcells.keys()
oldnames = oldcells.keys()
allnames = newnames + oldnames

for cell1 in allnames:
    for cell2 in allnames:
        if cell1 == cell2: continue
        dia_list1 = newcells[ cell1 ]
        dia_list2 = oldcells[ cell2 ] 
        N = max( len(dia_list1), len( dia_list2 ) )
        for i in range( N ):
            D1 = pd.Diagram( dia_list1 )
            D2 = pd.Diagram( dia_list2 )
            
            # TODO ~ 
            D1.write_diagram()
            D2.write_diagram()
            
            # also, must be able to take diagram as an argument
            dist = D1.compute_wasserstein_distance( D2 )

