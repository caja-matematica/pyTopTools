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
import matplotlib.pyplot as plt

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
        else:
            self.diagram_dim = dim

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
            fig = plot_diagram_scaled( self.data, scale=scale, fig=fig,
                                       inf_value=self.inf_value, **args )
        # if not scale:
        #     fig = pers.plot_diagram( fname, fig=fig, **args )
        # else:
           
        return fig


## utility functions


def find_unique( arr ):
    """
    arr : 1D array 
    """
    uniq = np.unique( arr )
    idx = [ np.where( arr==x )[0] for x in uniq ]
    return idx


def plot_diagram_scaled( diagram, fontsize=12, scale=None, color='b',
                         show_fig=True, fig=None, title=None, resize=False,
                         inf_value=None, marker_scale=1, show_legend=False, **args ):
    """
    This is the smae as plot_diagram(), except that each point on the
    diagram is scaled in relation to its number of occurences in the
    persistence diagram multiset. Thus, a birth-death pair will garner
    a larger marker if there are a realtively large number of unique
    generators with the same birth-death coordinates.
    
    persFile -- path to <perseus output>_*.txt, where * is the dimension.

    scale -- Factor to scale the birth/death times.

    marker_scale -- Factor by which to scale marker sizes.
    """
    if hasattr( diagram, '__array__' ):
        s = diagram
    else:
        if scale:
            # cast values as floats for division
            s = np.loadtxt( diagram, dtype=np.float, delimiter=' ' )
            s /= scale
        else:
            s = np.loadtxt( diagram, dtype=np.int, delimiter=' ' )
        
    try:
        births = s[:,0]
        deaths = s[:,1]
    except IndexError:
        # s is an (n,) array, so it must be reshaped in-place for
        # proper indexing
        print s.shape
        s.resize( ( s.shape[0], 1 ) )
        print s
        births = s[0] 
        deaths = s[1]

    # max death time
    if deaths.max() > 0:
        if inf_value==None:
            maxd = deaths.max()
        else: 
            maxd = inf_value
    else:
        maxd = births.max()
    print "Max death time ",  maxd

    # non-infinite gens
    normal_idx = np.where( deaths != -1 )[0]

    # find indices of unique birth-death coords.
    # [1] index just pulls off the indices
    uniq_idx = find_unique( deaths[normal_idx] )

    # add to an existing figure if necessary
    if not fig:
        fig = plt.figure( ) 
        fig.patch.set_alpha( 0.0 )
    try:
        ax = fig.gca()
    # in case we pass in the axes instance itself
    except AttributeError:
        ax = fig

    if len( normal_idx ) > 0:
        for u in uniq_idx:
            if resize:
                size = marker_scale*len( u )
            else:
                size = marker_scale
            ax.plot( births[ u ], deaths[ u ],
                     color+'o', ms=size, alpha=0.8 )
            # ax.plot( births[normal_idx], deaths[normal_idx], ,
            #          color+'o', )

    # create and plot the diagonal
    diag = [0, maxd+2]
    ax.plot(diag, diag, 'g-')

    # infinite gens
    inf_idx = np.where( deaths < 0 )[0]
    inf_vec = (maxd + 1) * np.ones( len( inf_idx ) )

    # plot the infinite generators
    ax.plot( births[inf_idx], inf_vec, 'ro', 
             label='Robust generators (num=' + str(len(inf_idx) ) + ')', 
             **args )
    # xticks = [ int( tk ) for tk in ax.get_xticks() ]
    # yticks = [ int( tk ) for tk in ax.get_yticks() ]
    ax.set_xticklabels( ax.get_xticks(), fontsize=fontsize )
    ax.set_yticklabels( ax.get_yticks(), fontsize=fontsize )
    ax.set_xlabel( 'birth', fontsize=fontsize )
    ax.set_ylabel( 'death', fontsize=fontsize )
    # fix the left x-axis boundary at 0
    ax.set_xlim( left=0, right=maxd+2 )
    ax.set_ylim( bottom=0, top=maxd+2 )

    # legend displaying number of robust/infinite generators
    if show_legend:
        ax.legend( loc=4 ) # 4 == lower right
    
    if title:
        ax.set_title( title, fontsize=16 )
    if show_fig:
        fig.show()
        
    print "Total number of persistence intervals", len( births ) 
    return fig


    
