"""
spatial.py

Author: Jesse Berwald

Opened: April 26, 2013

Spatial correlation analysis of the Ricker patch model.

Requirements: PySal (http://pythonhosted.org/PySAL/index.html)
"""
import pysal
import numpy as np

class Spatial:
    """
    Spatial correlation indicators, to replicate and extend results of
    Dakos.
    """
    def __init__( self, x, w=None ):
        """
        x : numpy array, either rectangular spatial array that will be
        raveled for pysal, or array that is already 1D. If x is 1D,
        then we must have that len(x)=k**2 in order to create the
        square weight matrix.

        w : if x is 1D, must provide w, a weights matrix indicating
        the spatial correlation between neighboring data points. Eg.,
        w_ij = 1 if grid element i is a neighbor of j, 0 otherwise.
        """
        
        if x.ndim > 1:
            # for ricker computations assume x is square
            n = x.shape[0]
            self.x = x.ravel()
               
            # create contiguous weight matrix ('rook' format, eg. NWSE
            # neighbors)
            self.w = pysal.lat2W( n, n )

        elif x.ndim == 1 and w is None:
            print "Must provide wieghts 'w'."
            return
        else:
            self.x = x
            n = int( np.sqrt( len( x ) ) )
            self.w = w

        # create the Moran object
        self.moran = pysal.Moran( self.x, self.w )
         
    def get_moran_I( self ):
        return self.moran.I


if __name__ == "__main__":

    # test using pysal's facilities
    f = pysal.open(pysal.examples.get_path("stl_hom.txt"))
    y = np.array(f.by_col['HR8893'])
    w = pysal.open(pysal.examples.get_path("stl.gal")).read()

    S = Spatial( y, w )
     
    # test square matrix input using random 5 x 5, make sure weights
    # created correctly
    #y = np.random.random( size=(5,5) )
    y = np.ones( (5,5) )
    y[0,0] = 0
    S2 = Spatial( y )

    
