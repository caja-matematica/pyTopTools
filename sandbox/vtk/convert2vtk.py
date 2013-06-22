import os, time
import numpy
from evtk.hl import imageToVTK
try:
    import chomp_betti as cb
except ImportError:
    raise

space = ' '

def vtk_coord( coord, dim=0 ):
    """
    Remove '(' and ')' from the points in the CUB file and return
    coord <dim>.
    """
    sc = coord.strip( '(' ).rstrip( ')\n' )
    point = sc.split( ',' )
    return space.join( point )+'\n'

def arr2vtk( coord ):
    """
    Convert array coord (text or .npy) to a vtk file.
    """
    point = space.join( map( str, coord ) )
    return point+'\n'

def cub2vtk( cubfile, fname, cubtype='chomp' ):
    """
    Write a VTK-readable file from a CUB file. We use 'Structured
    Grid' type for now. The model that we use is

    'octree_cylinder_in_cylinder.vtk'
    """
    if cubtype == 'chomp':
        # loop optimizations
        vc = vtk_coord
        # convert CUB format for each coordinate
        with open( cubfile, 'r' ) as fh:
            pts = map( vc, fh.readlines() )
    elif cubtype == 'arr':
        av = arr2vtk
        try:
            arr =numpy.loadtxt( cubfile, delimiter=' ', dtype='uint' )
            pts = map( av, arr )
        # catch .npy file format here (later)
        except IOError:
            raise
    # write VTK file 
    with open( fname, 'w' ) as fh:
        num_pts = len( pts )
        # write preamble
        # DIMENSIONS nx ny nz where nx=ny=nx in the cubicle case
        preamble = [ '# vtk DataFile Version 3.0\n',
                    cubfile+' vtk conversion\n',
                    'ASCII\n',
                    'DATASET STRUCTURED_GRID\n',
                    'DIMENSIONS 1 1 1\n',
                     'POINTS '+str( num_pts )+' double\n'
                     ]
        # print "preamble", preamble
        fh.writelines( preamble )
        # write coordinates
        map( fh.writelines, pts )
    return pts

def image2VTK( arr, vname=None ):
    """
    Save cubical or array 'image' data to VTK format.

    Requires: EVTK package
    """
    if type( arr ) == str:
        # set vtk prefix to same as array name
        vname = fname = arr
        # time shit
        tstart = time.time()
        print "Loading cubical complex coordinates..."
        if not fname.endswith( 'txt' ): fname+='.txt'
        # arr is now a numpy array
        arr = numpy.loadtxt( fname, dtype='int' )
        nx, ny = arr.shape
        print "Loaded! shape=(", nx, ny, ")"
        print "Elapsed time:", time.time() - tstart
    else:
        if vname == None:
            raise ValueError, "Must provide a name for VTK file if passing in numpy array!"
 
    # shift the origin to (0,0,0)
    for i in range( ny ):
        arr[:,i] = arr[:,i] - arr[:,i].min()
    dx, dy, dz = [ int( arr[:,i].max() ) for i in range( ny ) ]
    # total number of cells in structured grid
    ncells = (dx+1) * (dy+1) * (dz+1)
    cells = numpy.zeros( ncells, dtype='uint' ).reshape( (dx+1,dy+1,dz+1), order='C' )

    print "Converting cubical structure to VTK grid structure..."
    print "Size of cell complex", cells.shape
    tstart = time.time()
    w = numpy.where
    for row in arr:
        cells[ row[0], row[1], row[2] ] = 1
    imageToVTK("./"+vname, cellData = {"cell" : cells} )
    print "Elapsed time:", time.time() - tstart

def cub2arr( cubfile, dim=3 ):
    """
    Convert a CUB file to a numpy array and write to text file, one
    coordinate per line.
    """
    print "Loading cubical complex..."
    return cb.pix2array( cubfile+'.cub', dim=dim )
    

if __name__ == "__main__":

    cubname = 'new_stack_1000f'  #'cub_large.txt'
    #fname = 'new_stack
    # first, convert CUB file to txt array'
    arr = cub2arr( cubname )
    print "Saving text file..."
    numpy.savetxt ( cubname+'.txt', arr, fmt='%d' )
    image2VTK( cubname )
    
