"""
perisstence_diagrams.py

author: Jesse Berwald

opened: September 20, 2013

Container for operations on persistence diagrams. Eg.,

- bottleneck and wasserstein distances
- diagram plotting

"""
import numpy as np
import subprocess as sp
import tempfile

# local stuff
from pyTopTools import perseus_wrap as pers


class Diagram( object ):
    """
    Base class. 

    input: Perseus diagram output at a given
    dimension. filename_'dim'.txt
    """
    def __init__( self, diagram, dim=None, inf=None ):
        """
        diagram : full path to diagram file on disk, or array
        containing diagram information.

        dim : [optional] We can extract this from the filename if
        necessary.
        """
        # diagram is an array, so save a temp file for distance
        # functions to read
        if hasattr( diagram, '__array__' ):
            fh = tempfile.NamedTemporaryFile(delete=False)
            np.savetxt( fh, diagram )
            self.data = diagram
            self.diagram = fh.name
        else:
            self.diagram = diagram
            # in case we need the data for other applications
            self.data = np.loadtxt( diagram )
        
        # the 'infinite' value usually comes from the maximum number
        # of steps in growing epsilon balls (say). If we don't have
        # this, then set inf to be the max death time +1
        if inf == None:
            self.inf_value = max( self.data[:,1] ) + 1

        if dim is None:
            idx = self.diagram.rfind( '.' )
            self.diagram_dim = int( self.diagram[idx-1] )

    def __repr__( self ):
        s = "Persistence diagram with "+ str( len( self.data ) )+\
            " generators in dimension "+ str( self.diagram_dim )
        return s

    def compute_bottleneck_distance( self, other ):
        """
        Compute the bottleneck distance between this the persistence
        diagram for this data and a diagram for 'other'.

        Diagram for this Window object must have been computed already.

        this_dia : path to this window's persistence diagram data
        (should be stored in self.persout attribute)

        other : path to persdia file

        engine refers to either Miro's code ('c') or Kelly's
        ('python'). Not implemented.
        """
        this_dia = self.diagram
        # grab the filename of the diagram if necessary
        if hasattr( other, 'diagram' ):
            other = other.diagram 
        try:
            dist = sp.check_output( ["bottleneck", this_dia, other] )
        except:
            print "subprocess returned an error!"
        return float( dist )

    def compute_wasserstein_distance( self, other ):
        """Compute the Wasserstein distance between self.persdia and other.

        Diagram for this Window object must have been computed
        already.

        this_dia : path to this window's persistence diagram data
        (should be stored in self.persout attribute)

        other : path to persdia file
        """
        this_dia = self.diagram
        # grab the filename of the diagram if necessary
        if hasattr( other, 'diagram' ):
            other = other.diagram 
        try:
            dist = sp.check_output( ["wasserstein", this_dia, other] )
        except:
            print "subprocess returned an error!"
        return float( dist )
        
    def draw_diagram( self, fig=None, scale=1, dim=None, **args ):
        """fname : full path to persistence diagram file. 
        
        fig : Figure object, in case we want to plot diagram on top of
        one another (for some reason...)

        scale : This scales (birth,death) pairs. Written for Morse
        function version of persistence. But useful in that the
        plot_diagram_scaled() function also plots markers whose size
        is determined by the number of generators at that (b,d)
        coordinate.

        """
        if scale:
            fig = pers.plot_diagram_scaled( self.diagram, scale=scale, fig=fig,
                                            inf_value=self.inf_value, **args )
        # if not scale:
        #     fig = pers.plot_diagram( fname, fig=fig, **args )
        # else:
           
        return fig
