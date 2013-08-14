#!/usr/bin/python
"""
ricker_patch_model.py

Author: Jesse Berwald

Opened: Match 26, 2013
"""
import os
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from pyTopTools import process_perseus as proc_pers
from pyTopTools import persdia_sub_super as pd
from pyTopTools import chomp_tools as chomp
from pyTopTools.ricker import spatial

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
        self.iter_number = 0
        self._square = square
        self._verbose = verbose
        self._IC = IC
        # change flag if creating dispersal matrix. 
        self._do_dispersal = False 

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

    def threshold( self, thresh=-1 ):
        """
        Threshold the population matrix at 'thresh'.

        If thresh == -1, use ln( fitness ).
        """
        if thresh == -1:
            thresh = np.log( self.fitness )
        # find all entries above the threshold
        ma = np.ma.masked_less( self.population, thresh )
        # return only the T/F mask
        return ma.mask

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
        self.iter_number += 1 # keeps track for file names
            
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

    def binarize( self, thresh=-1 ):
        """
        Threshold to get a binary array, with 1's where
        self.population >= thresh. See threshold() for usage of
        thresh.
        """
        mask = self.threshold( thresh )
        z = np.ones_like( self.population, dtype=int )
        z[mask] = 0
        return z

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
        binary = self.binarize()
        cub = chomp.array2cub( binary )
        chomp.write_cubical_file( cub, chomp_in )
        out = chomp.run_chomp( chomp_in )
        betti = chomp.extract_betti_string( out )
        return betti

    def compute_persistence( self, persname, scale, shift_mean=False ):
        """
        Compute sub- and superlevel set cubical topological
        persistence on population matrix.

        shift_mean : this will shift many values below zero, which we
        probably do not want.
        """
        if shift_mean:
            pop = self.population - self.population.mean()
        else:
            pop = self.population
        proc_pers.perseus_sub_super( pop, persname, scale )

    def compute_moran( self ):
        """
        Computes Moran's coefficient for (local) spatial correlation.
        """
        S = spatial.Spatial( np.asarray(self.population) )
        return S.get_moran_I()
        

    def draw_persistence_diagram( self, fname, **args ):
        """
        """
        if not os.path.isfile( fname ):
            print "Persistence diagram file ", fname, " does not exist! "+\
                "Please create it using self.compute_persistence()."
            return
        else:
            fig = pd.persdia_sub_super( fname, **args )
            return fig

    def draw_matrix( self, thresh=None, val=-1 ):
        """
        Use imshow() to to display the population matrix.
        """
        if thresh:
            m = self.threshold( val )
        else:
            m = self.population
        fig = plt.figure()
        ax = fig.gca()
        ax.imshow( m, interpolation='nearest' )
        return fig

    def draw_matrix_3d( self ):
        """
        """
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
    
    # ========= END RICKER CLASS ===========

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

def run( IC, tfinal, fit, disp, spin_up=0, ndim=1,
         square=True, save_pop=False, save_abundance=False,
         betti_name=None, pers_name=None, pers_scale=1, compute_moran=False ):
    """
    betti_name : prefix for chomp .cub and betti files

    pers_name : path and name of persistence files for Perseus.

    pers_scale : optional scaling of population matrix in persistence
    computation. Necessary for Perseus's integer requirements.
    
    """
    R = Ricker( IC, fitness=fit, dispersal=disp, ndim=ndim, square=square )

    # containers for storage if necessary
    if save_pop:
        pop = [ R.population ]
        binary = [ R.binarize() ]
    if save_abundance:
        abundance = [ R.abundance() ]
    if betti_name:
        betti = [ R.compute_betti( betti_name + str(R.iter_number) ) ]
    if compute_moran:
        moran = [ R.compute_moran() ]

    # Now iterate system
    if spin_up > 0:
        for i in range( spin_up ):
            R.population_map()            
    for i in range( spin_up, tfinal ):
        R.population_map()
        if save_pop:
            pop.append( R.get_population_matrix() )
            binary.append( R.binarize() )
        if save_abundance:
            abundance.append( R.abundance() )
        if betti_name:
            betti.append( R.compute_betti( betti_name + str(R.iter_number) ) )
        if pers_name:
            R.compute_persistence( pers_name+str(R.iter_number),
                                   pers_scale )
        if compute_moran:
            moran.append( R.compute_moran() )

    output = { 'ricker': R }

    if save_pop:
        output['population'] = pop
        output['binary'] = binary
    if save_abundance:
        output['abundance'] = abundance
    if betti_name:
        b = np.asarray( betti )

        print "b", b
        
        betti = b[:,:2]  # just keep 0 and 1 homology levels
        output['betti'] = betti
    if compute_moran:
        output['moran'] = moran
    return output
    

if __name__ == "__main__":

    IC='./debug/test3x3.txt'
    disp = 0.0
    fit = 22
    ndim = 20
    tf = 20
    pers_scale = 1000
    betti_path ='debug/rickerBETTI_TEST' #betti = 'debug/bettiTEST'
    
    #ic = init_population( IC )
    ic = 10*np.random.random( (ndim,ndim) )

    out = run( ic, tf, fit, disp, spin_up=10, betti_name=betti_path,
               save_abundance=True, save_pop=True, pers_scale=pers_scale,
               pers_name='debug/rickerPERS_TEST', compute_moran=True )

    #persname = 'debug/persTEST'
    #proc_pers.perseus_sub_super( R.population, persname, 1000 )

    
