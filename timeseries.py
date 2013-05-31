import numpy as np
import subprocess as sp
import matplotlib.pyplot as plt
import npy2perseus as n2p
from pyTopTools import perseus_wrap as pers
from pyTopTools import bottleneck_distance as BD
import tempfile

"""
Module for handling time series and ndarrays for persistence analysis.

Author: Jesse Berwald

Opened: March 10, 2013
"""

class Timeseries( object ):
    """
    Basic container for data. Methods for reading, writing, analyzing.

    data : numpy array or path to file (.txt or .npy)

    data_type : cubical ('cub'), sparse cubical ('scub'), or
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

        self.dtype = data_type


    def __repr__( self ):
        s = "Timeseries with " + str( len(self.data) ) +\
            " points. min = "+str(self.data.min())+", max = "+\
            str( self.data.max() )
        return s

    def mean( self ):
        return self.data.mean()

    def var( self ):
        return self.data.var()

    def std( self ):
        return self.data.std()

    def convert2perseus( self, persname, **kwargs ):
        """
        Convert an array, or text or numpy file to perseus format.

        data -- numpy array

        Writes a 1D point cloud to Persues-readable Vietoris-Rips format.
        """
        space = " "
        fargs = {'radius_scaling' : 1,
                 'stepsize' : 0.02,
                 'nsteps' : 10,
                 'bradius' : 0.0,
                 'embed_dim' : 1
                 }
        fargs.update( kwargs )

        # fix up end of filename
        fargs['persname'] = persname
        if not fargs['persname'].endswith( 'txt' ):
            fargs['persname'] += '.txt'

        # This gets appended to the end of every line
        br = str( fargs['bradius'] )

        with open( fargs['persname'], 'w' ) as fh:
            # ambient dimension
            fh.write( str( fargs['embed_dim'] )+'\n' )

            # initial threshold, step size, num steps, dimension cap==ny
            params = [ str( fargs['radius_scaling'] ),
                       str( fargs['stepsize'] ),
                       str( fargs['nsteps'] )
                       ]
            params = space.join( params )
            params += '\n'
            fh.write( params )

            # now write the timeseries and birth radii to disk
            for obs in self.data:
                try:
                    r = [ str( x ) for x in obs ]
                except TypeError:
                    # quicker than if/else 
                    r = [ str( obs ) ]
                r += [ br ] # add the birth radius
                r = space.join( r )
                r += '\n'
                fh.write( r )
        print "wrote file to ", fargs['persname']
        out = { 'filename' : fargs['persname'],
                'data' : self.data }
        return out

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

class Window( Timeseries ):

    def __init__( self, data, tmin=0, tmax=-1, data_type='vr' ):

        Timeseries.__init__( self, data, data_type )
        self.tmin = tmin
        self.tmax = tmax
        if not tmin:
            tmin = 0
        if not tmax:
            tmax = len( data )
        self.data = self.data[ tmin:tmax ]
        self.persdia = None

    def compute_persistence( self, fname, output=None, dtype='brips', debug=False ):
        """
        Compute persistence topology of data cloud using Perseus on
        the VR complex in fname.
        """
        # prefix for the persistence diagram files
        if not output:
            if fname.endswith( '.txt' ):
                output = fname[:-4]
        pers.perseus( fname, output, dtype, debug=debug )

        # just load 0-dim diagram for now since we're just working with time series
        self.persdia = np.loadtxt( output + '_0.txt' )

    def compute_bottleneck_distance( self, other, this_dia=None,
                                     return_match=False, engine='py' ):
        """
        Compute the bottleneck distance between this the persistence
        diagram for this data and a diagram for 'other'.
        """
        # the persistence diagram for this window has not already been computed.
        if self.persdia is None:
            try:
                self.persdia = np.loadtxt( this_dia )
            except IOError:
                return
            self.compute_persistence( self.persdia )
        if engine =='py':
            dist, matching = BD.bdist( self.persdia, other )
            if not return_match:
                return dist
            else:
                return dist, matching
        # insert call to Miro's C++ bottleneck code
        elif engine=='c':
            try:
                p = sp.check_output( ["bottleneck", this_dia, other] )
            except:
                print "subprocess returned an error!"
            return int( p )

    def compute_wasserstein_distance( self, this_dia, other ):
        """
        Compute the Wasserstein distance between self.persdia and other.

        this_dia : path to this window's persistence diagram data

        other : path to persdia file
        """
        try:
            p = sp.check_output( ["wasserstein", this_dia, other] )
        except:
            print "subprocess returned an error!"
        return int( p )
        
    def draw_diagram( self, fname=None, dim=0, fig=None, scale=False, **args ):
        """
        """
        if not fname:
            fname = self.persdia
            fname = fname + '_' + str( dim ) + '.txt'
        if not scale:
            fig = pers.plot_diagram( fname, fig=fig, **args )
        else:
            fig = pers.plot_diagram_scaled( fname, fig=fig, **args )
        return fig
    
    
