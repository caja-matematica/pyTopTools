#!/usr/bin/python
"""
ricker_patch_model.py

Author: Jesse Berwald

Opened: Match 26, 2013
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from pyTopTools import process_perseus as proc_pers
from pyTopTools import persdia_sub_super as pd
from pyTopTools import chomp_tools as chomp

class Ricker( object ):
    """
    Basic object to run the deterministic Ricker map on a square
    population matrix.
    """
    def __init__( self, IC, fitness=5, dispersal=0, ndim=1, square=False,
                  verbose=False ):
        """
        IC -- initial condition. A scale, list, or matrix, depending
        on dimension of system (see ndim)

        fitness -- The fitness on each patch. If scalar, each patch
        has the same fitness. If numpy matrix of fitness values (same
        dimension as ndim), different patches can have different
        fitness values.

        dispersal -- Dispersal parameter. As with fitness, can be either
        scalar or numpy matrix.

        ndim -- dimension of the system.

        verbose -- output debugging and other chatter. (default=False)
        """
        self.dispersal = dispersal # possibly a vector of matrix of different values
        self.fitness = fitness # ditto
        self.ndim = ndim
        self._square = square
        self._verbose = verbose
        self._IC = IC
        self._do_dispersal = False # change flag if creating dispersal matrix

        if square:
            self.sysdim = ndim * ndim
        else:
            self.sysdim = ndim

        # now create patch population container
        if not hasattr( IC, '__iter__' ):
            # single initial population value will fill entire pop. matrix
            if self._square:
                self.population = np.empty( (ndim,ndim) )
                self.population[:] = IC
            # for single patch case
            else:
                self.population = IC
        else:
            # not much error catching here.
            self.population = np.matrix( IC )

        if self._square:
            # create matrix of off-diagonals for dispersal computation
            if not hasattr( self.dispersal, '__array__' ):
                d = np.ones( self.ndim - 1 ) # clip one entry to account for
                                             # off-diag position
                self.D = np.diagflat( d, 1 ) + np.diagflat( d, -1 )
                self.D = np.matrix( self.D ) # for matrix multiplication
                self._do_dispersal = True
            else:
                print "Multiple values for dispersal not implemented yet!"

    def __repr__( self ):
        s = "Ricker object with " + str(self.sysdim) +" patches, fitness "+\
            str(self.fitness)+" and dispersal "+str(self.dispersal)+"."
        return s

    def threshold( self, thresh ):
        """
        Threshold the population matrix at 'thresh'.
        """
        ma = np.ma.masked_less( self.population, thresh )

    def abundance( self ):
        try:
            return self.population.sum()
        # in case abundance is called for a simulation with a single
        # patch.
        except AttributeError:
            return self.population

    def get_population_matrix( self ):
        return self.population

    def population_map( self, return_values=False ):
        """
        Perform one iteration of the Ricker map. Update
        self.population. To access, see get_populaton_matrix().
        """
        # update population on each patch. Call np.multiply for
        # elementwise multiplication
        P = self.fitness * np.multiply( self.population, np.exp( -self.population ) )
        if self._do_dispersal:
            self.population = ( 1 - self.dispersal ) * P + \
                  ( self.dispersal / 4. ) * ( P * self.D ) + \
                  ( self.dispersal / 4. ) * ( self.D * P )
        else:
            self.population = P
            
    def save_pop( self, output_name ):
        """
        output_name : write data to 'output_name' (extension
        '_pop.npy' and '_abun.npy' appended to filename).
        """
        if self._verbose:
            print "Saving files..."
        # stack the array in time
        pop = np.asarray( self.population )
        np.save( output_name + '_pop.npy', self.population )

        # if tracking abundance data, save it!
        abun = np.asarray( self.abundance )
        np.save( output_name + '_abun.npy', abun )

        if self.verbose:
            print "Done!"

    def compute_betti( self, output ):
        """
        Calls chomp_tools.run_chomp() on output.

        output -- where to write the chomp-readable file.
        """
        if not output.endswith( 'cub' ):
            chomp_in = output + '.cub'
            chomp_out = output + '.cbetti'
        else:
            chomp_out[:-4] + '.cbetti'
        chomp.array2chomp( self.population, chomp_in )
        out = chomp.run_chomp( chomp_in )
        betti = chomp.extract_betti_string( out )
        return betti

    def compute_persistence( self, persname, scale ):
        """
        Compute sub- and superlevel set cubical topological
        persistence on population matrix.
        """
        proc_pers.perseus_sub_super( self.population, persname, scale )

    def draw_persistence_diagram( self, fname ):
        """
        """
        if not os.path.isfile( fname ):
            print "Persistence diagram file ", fname, " does not exist! "+\
                "Please create it using self.compute_persistence()."
            return
        else:
            fig = pd.persdia_sub_super( fname )
            return fig

    def draw_matrix( self ):
        """
        Use imshow() to to display the population matrix.
        """
        fig = plt.figure()
        ax = fig.gca()
        ax.imshow( self.population )
        return fig
    
    # ========= END RICKER ===========

def init_population( IC ):
    """
    If filename, read an IC population matrix from
    disk. Otherwise, just append to self.population.
    """
    try:
        ic = np.load( IC ) 
    # seems that IC was not a path to a .npy file, try loadtxt()
    except IOError:
        ic =np.loadtxt( IC ) 
    return ic

def run( IC, tfinal, fit, disp, square=True, save_abundance=False,
         betti_name=None ):
    """
    betti_name -- prefix for chomp .cub and betti files
    """
    R = Ricker( ic, fitness=fit, dispersal=disp, ndim=ndim, square=square )

    pop = []

    if save_abundance:
        abundance = []
    if betti_name:
        betti = []
    for i in range( tfinal ):
        R.population_map()
        pop.append( R.get_population_matrix() )
        if save_abundance:
            abundance.append( R.abundance() )
        if betti_name:        
            betti.append( R.compute_betti( betti_name ) )

    output = [ R, pop ]
    if save_abundance:
        output.append( abundance )
    if betti_name:
        output.append( betti )
    return output
    

if __name__ == "__main__":

    IC='./debug/test3x3.txt'
    disp = 0.0
    fit = 5
    ndim = 3
    tf = 10
    betti = 'debug/bettiTEST'
    
    ic = init_population( IC )
    #ic = np.random.random( (ndim,ndim) )

    out = run( ic, tf, fit, disp, betti_name=betti, save_abundance=True )

    persname = 'debug/persTEST'
    #proc_pers.perseus_sub_super( R.population, persname, 1000 )

    
