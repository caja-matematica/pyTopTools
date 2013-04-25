import numpy as np
import subprocess as sp
import matplotlib.pyplot as plt
import npy2perseus as n2p

"""
Module for handling time series and ndarrays for persistence analysis.

Author: Jesse Berwald

Opened: March 10, 2013
"""

class Timeseries( object ):
    """
    Basic container for data. Methods for reading, writing, analyzing.

    data : numpy array or path to file (.txt or .npy)

    data_type : cubical ('cub'), sparse cubical ('scub'),
    vietoris-rips point cloud ('vr'). 
    """
    def __init__( self, data, data_type ):

        if hasattr( data, '__array__' ):
            self.data = data
        else:
            self.fname = data
            try:
                self.data = np.load( self.fname )
            except IOError:
                self.data = np.fromfile( self.fname, sep='\n' )

                # self.dtype = data_type


    def convert2perseus( self, persname, **kwargs ):
        """
        Convert an array, or text or numpy file to perseus format.
        """
      #  if dtype == 'timeseries':
        out = n2p.write_timeseries( self.data, persname, **kwargs )
        # elif dtype == 'scub':
        #     out = n2p.write_scubtop( self.data, persname, **kwargs )
        # elif dtype == 'cub':
        #     out = n2p.write_cubtop( self.data, persname, **kwargs )
        # elif dtype == 'vr':
        #     out = n2p.write_vr( self.data, persname, **kwargs )
        # else:
        #     print "Unknown data type!"
       

    def draw_data( self, **kwargs ):
        """
        Either plot the data as a time series, a scatter plot (no
        implemented) or use imshow for a matrix.
        """
        if self.data.ndim == 1:
            fig = plt.figure()
            ax = fig.gca()
            ax.plot( self.data, **kwargs )
            return fig
        elif self.dtype == 'scub' or self.dtype == 'cub':
            fig = plt.figure()
            ax = fig.gca()
            ax.imshow( self.data )
            return fig

    
