"""
Module to convert numpy array to various Persues formats.
"""

import numpy
import matplotlib.pylab as plt
import re
import os
import pickle as pkl
from scipy.spatial import distance

space = ' '  # for strings

def write_file ( fname, output ):
    """ 
        ----DENSE FORMAT---
        - Write contents of single frame to output file
            e.g. new cell name on simplex:
            '/data/jberwald/wyss/data/Cells_Jesse/New/frames/new_110125/new_110125-concatenated-ASCII_324.npy'
            e.g. old cell name on simplex:
            '/data/jberwald/wyss/data/Cells_Jesse/Old/frames/old_100125/old_100125-concatenated-ASCII_1525.npy'
        Perseus input requires following format:
        #Dimension of cubical grid
        #Number of cubes in first dimension
        #Number of cubes in second dimension
        #Number of cubes in   ...     dimension
        #Birth time of cube anchored at (0,0,...,0)
        #Birth time of cube anchored at (0,0,...,1)
        #Note: Output is file name
    """
    #Load npy file
    data = numpy.load(fname)
    size = data.shape#(203,198)
    dim = len ( size )
    #Open file
    text = open (output + ".txt", "w")
    #Write number of dimensions
    text . write (str(dim) + "\n")
    #Write resolutions
    for i in xrange(len(size)):
        text . write (str(size[i]) + "\n")
    #Write strictly 2D file
    for i in xrange(size[0]):#203
        for j in xrange(size[1]):#198
            if int(data[i][j]) == 0:
                text . write (str(-1) + "\n")
            else:
                text . write (str(int(data[i][j])) + "\n")
    text . close()

def write_timeseries( arr, output, scale=1, block=None, ptype='cub', **kwargs ):
    """
    Write blocks of timeseries data to disk in specified format.

    arr : numpy 1d array.

    output : Path to perseus-readable file.

    scale : (1, \infty), default = 1. If scale != 1, values of arr are
    rounded up to ceiling.

    block : tuple (start,end) of indices to slice and write to disk (instead of complete arr).
    """
    if block:
        data = numpy.copy( arr[ block[0]:block[1] ] )
    else:
        data = numpy.copy( arr )
        
    min_val = data.min()
    if min_val < 1:
        data += abs( min_val ) + 1

    # For *cub formats it's necessary to rescale if the data are
    # floats, and then take the ceil (or floor) by converting to ints.
    if scale != 1:
        data *= scale
        data = numpy.asarray(data, dtype=int )

    if ptype == 'cub':
        with open( output, 'w' ) as fh:
            fh.write( '1\n' ) # assume 1D data
            fh.write( str( data.shape[0] ) + '\n' )
            for x in data:
                fh.write( str( x ) + '\n' )
    elif ptype == 'vr':
        write_vr( data, output=output, **kwargs )

def write_cubtop( arr, output, ndim=2, scale=1, dtype=None ):
    """
    Write an array of values to dense cubical toplex Perseus format.

    Note: Only implemented for dims 1 and 2.
    """
    min_val = arr.min()
    if min_val < 1:
        arr += abs( min_val ) + 1

    # For *cub formats it 's necessary to rescale and take the ceil
    # (or floor)
    if scale != 1:
        arr *= scale
    if dtype:
        arr = numpy.asarray( arr, dtype=dtype )    
        
    with open( output, 'w' ) as fh:
        fh.write( str( arr.ndim )+'\n' )
        # row and column dims for Perseus memory alloc
        for d in range( arr.ndim ):
            fh.write( str( arr.shape[d] ) + '\n' )

        if arr.ndim == 1:
            for x in arr:
                fh.write( str( x ) + '\n' )
        else:
            # write every value to disk
            for i in xrange( arr.shape[0] ): 
                for j in xrange( arr.shape[1] ):
                    fh.write( str( int( arr[i,j] ) ) + '\n' )

    return output

def write_scubtop( arr, output, ndim=2 ):
    """
    Write an array to *sparse* Perseus format.
    
    arr -- numpy ndarray

    output -- name (full path) of output file

    NOTE -- Only implemented for arr.ndim == 2.
    """
    with open( output, 'w' ) as fh:
        fh.write( str( arr.ndim )+'\n' )
        for i in xrange( arr.shape[0] ): 
            for j in xrange( arr.shape[1] ):
                if int( arr[i,j] ) != 0:
                    pos = str(i) + ' ' + str(j) + ' '
                    fh.write (pos + str(int( arr[i,j] )) + "\n")
    return output


def write_vr( data, **kwargs ):
    """
    data -- numpy array

    Writes a point cloud to a Persues-readable Vietoris-Rips format.
    """
    fargs = {'output' : None,
             'radius_scaling' : 1,
             'stepsize' : 0.2,
             'nsteps' : 10,
             'bradius' : 0.1,
             'tstepsize' : 1,
             'embed_dim' : 1
             }
    fargs.update( kwargs )

    if not fargs['output'].endswith( 'txt' ):
        fargs['output'] += '.txt'

    # This gets appended to the end of every line
    br = str( fargs['bradius'] )

    with open( fargs['output'], 'w' ) as fh:
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
        for obs in data:
            try:
                r = [ str( x ) for x in obs ]
            except TypeError:
                # quicker than if/else 
                r = [ str( obs ) ]
            r += [ br ] # add the birth radius
            r = space.join( r )
            r += '\n'
            fh.write( r )
    print "wrote file to ", fargs['output']
    out = { 'filename' : fargs['output'],
            'data' : data }
    return out


def write_distance_mat( data, output=None, g=0.1, stepsize=0.2,
                        nsteps=10, dcap=None, metric='euclidean' ):
    """
    For use with

    perseus distmat <path to distance matrix file> <output string>
    """
    # load from file if data is not an array of points
    if not hasattr( data, '__index__' ):
        # .npy or .txt
        try:
            data = numpy.load( data )
        except IOError:
            data = numpy.loadtxt( data )
    dist = distance.pdist( data, metric=metric )
    if output:
        dist = distance.squareform( dist )
        space = " "
        # num points x dimension
        nx, ny = data.shape
        if not dcap:
            dcap = ny
        if not output.endswith( 'txt' ):
            output += '.txt'
        with open( output, 'w' ) as fh:
            # nx x nx distance matrix
            fh.write( str( nx )+'\n' )
            # initial threshold, step size, num steps, dimension cap==ny
            params = [ str( g ), str( stepsize ), str( nsteps ), str( dcap ) ]
            params = space.join( params )
            params += '\n'
            fh.write( params )
            for row in dist:
                r = [ str( x ) for x in row ]
                r = space.join( r )
                r += '\n'
                fh.write( r )
        return dist
    else:
        return distance.squareform( dist )
