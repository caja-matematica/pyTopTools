import os, time
import numpy
from evtk.hl import imageToVTK
try:
    import chomp_tools as ct
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

def data2vtk( datafile, fname, dtype='arr', metadata='' ):
    """Write a VTK-readable file from a CUB file. We use 'Structured
    Grid' type for now. The model that we use is

    'octree_cylinder_in_cylinder.vtk'

    datafile : path to numpy array on disk (txt or npy); or path to
    CUB file. Also possible to pass in just the array itself.

    fname : name for the output (.vti, vtk, etc)

    dtype : 'arr' (default) or 'chomp'. eg., numpy array or cubical complex.

    """
    if dtype == 'chomp':
        # loop optimizations
        vc = vtk_coord
        # convert CUB format for each coordinate
        with open( datafile, 'r' ) as fh:
            pts = map( vc, fh.readlines() )
    elif dtype == 'arr':
        av = arr2vtk
        try:
            arr =numpy.loadtxt( datafile, delimiter=' ', dtype='uint' )
            pts = map( av, arr )
        # catch .npy file format here (later)
        except IOError:
            raise
        # we can just pass in the array itself
        except ValueError:
            arr = datafile
            pts = map( av, arr )
            
            print "PTS", pts[:10]

    # write VTK file 
    with open( fname, 'w' ) as fh:
        num_pts = len( pts )
        # write preamble
        # DIMENSIONS nx ny nz where nx=ny=nx in the cubicle case
        preamble = [ '# vtk DataFile Version 3.0\n',
                    metadata+' vtk conversion\n',
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
    Save array 'image' data to VTK format.

    Requires: EVTK package
    """
    if type( arr ) == str:
        # set vtk prefix to same as array name
        vname = fname = arr
        # time shit
        tstart = time.time()
        print "Loading array..."
        if not fname.endswith( 'txt' ): fname+='.txt'
        # arr is now a numpy array
        arr = numpy.loadtxt( fname, dtype='int' )
        print "Time to load file:", time.time() - tstart
    else:
        if vname == None:
            raise ValueError, "Must provide a name for VTK file if passing in numpy array!"
    nx, ny = arr.shape
    print "Loaded! shape=(", nx, ny, ")"
    print "shape", arr.shape
    # shift the origin to (0,0,0)
    # for i in range( ny ):
    #     arr[:,i] = arr[:,i] - arr[:,i].min()

    # arr records corners of cubical grid (structured grid) in
    # 3-space; the extent of the grid is found by the max value
    # in each dimension.
    dx, dy, dz = [ int( arr[:,i].max() ) for i in range( ny ) ]
    print "dx, dy, dz:", dx, dy, dz
    # total number of cells in structured grid
    ncells = (dx+1) * (dy+1) * (dz+1)
    cells = numpy.zeros( ncells, dtype='uint' ).reshape( (dx+1,dy+1,dz+1), order='C' )
    #cells = numpy.ones( ncells, dtype='uint' ).reshape( (dx+1,dy+1,dz+1), order='C' )

    print "Converting cubical structure to VTK grid structure..."
    print "Size of cell complex", cells.shape
    tstart = time.time()
    cells[ arr[:,0], arr[:,1], arr[:,2] ] = 1
    print cells
    # for row in arr:
    #     print row
    #     cells[ row[0], row[1], row[2] ] = 1
    imageToVTK( vname, cellData = {"cells" : cells} )
    print "Wrote", vname, "to file."
    print "Time to write vtk file:", time.time() - tstart

def cub2arr( cubfile, dim=3 ):
    """
    Convert a CUB file to a numpy array, one coordinate per line.
    """
    print "Loading cubical complex..."
    return ct.cub2array( cubfile+'.cub', dim=dim )
    

if __name__ == "__main__":

    start = time.time()
    # no .cub extension (added later)
    cubname = 'new_stack_1000f.txt'
    cubdir = './sandbox/vtktest/' #new_stack_1000f'  #'cub_large.txt'
    fname = cubdir + cubname

    #first, convert CUB file to txt array
    arr = cub2arr( fname )
   
    savedir = './sandbox/vtktest/'
    #numpy.savetxt ( savedir + cubname+'.txt', arr, fmt='%d' )
  #  image2VTK( arr, vname=savedir + cubname )
    
    print "Elapsed time: ", time.time() - start
