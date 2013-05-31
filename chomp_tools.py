import subprocess, os
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import numpy# as np
import re
import cPickle as pkl
import tempfile
try:
    import pp
except ImportError:
    print "No parallel python install"
    print "We'll try to continue..."
import time
from collections import deque
#from scipy.stats import linregress
# Python image libary
try:
    from PIL import Image  
except ImportError:
    print "Python Image Library not installed"
    raise

slash = '/'

def natural_key(string_):
    """
    Use with frames.sort(key=natural_key)
    """
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string_)]

def run_chomp( fname ): #, savename ):
    """
    Call chomp to compute the betti numbers of the image in file fname.

    See http://chomp.rutgers.edu
    """
    try:      
        p = subprocess.check_output( ["chomp-rutgers", fname] )
    except:
        print "subprocess returned with command: ", cmd
    return p

def array2cub( arr ):
    """
    arr : binarized (thresholded) array.
    
    Returns a list of indices suitable for writing a cubical complex.
    A cubical complex file consists of a list of many lines of the
    form ( n1, n2, ..., nd ) where d is the dimension of the complex
    and the ni's are the coordinates of the individual cube.

    Example: (0, 0, 1, 0) (10, 13, 2, 3) ... two cubes of a four
    dimensional cubical complex

    Note: This only works for 2D complexes
    """
    # sometimes np.matrix causes problems or unintuitive dimension
    # issues
    arr = numpy.asarray( arr )
    # find appropriate cubical corners
    w = numpy.where( arr==1 )
    return numpy.array( zip( w[0], w[1] ), dtype=int )

def write_cubical_file( cub, fname ):
    """
    cub : coordinates of cubical corners.

    fname : full path to output file. 
    """
    if not fname.endswith( '.cub' ):
        fname += '.cub'

    with open( fname, 'wb' ) as fh:
        for row in cub:
            x, y = str(row[0]), str(row[1])
            line = '(' + x + ',' + y + ')' + '\n'
            fh.write( line )
                   
def extract_betti_string( chomp_out ):
    """
    chomp_out -- string output from chomp-rutgers
    """
    out = chomp_out
    sout = out.split( '\n' )
    for line in sout:
        if line.startswith( 'Betti' ):
            # keep only the numbers
            # index=2 excludes the text
            betti_numbers = line.strip().split()[2:]
    return [ int( b ) for b in betti_numbers ]
                    
def extract_betti( fname, betti_file=None ):
    """
    Read the betti numbers from the 'fname' containing the output from
    chomp. Save raw Betti numbers
    """
    # open and read chomp-produced data file
    with open( fname, 'r' ) as fh:
        lines = fh.readlines()
    # grab the line with the Betti numbers
    if len( lines ) > 1:
        for line in lines:
            if line.startswith( 'Betti' ):
                # keep only the numbers
                betti_numbers = line.strip().split()[2:]
    # chomp output single line of betti number (Shaun Harker's "other"
    # chomp version...)
    else:
        betti_numbers = lines[0].strip().split()
    max_dim = len( betti_numbers )
    # open betti file and save generators. write the same format for
    # both chomp outputs to keep things consistent
    if not betti_file:
        betti_file = fname[:-3] + 'betti'
    with open( betti_file, 'w' ) as fh:
        for i, b in enumerate( betti_numbers ):
            line = str(i) + ' ' + betti_numbers[i] +'\n'
            fh.write( line )

def read_betti_dir( fdir, suffix='.hom' ):
    """
    Read all .betti files in a directory and organize them for analysis.
    """
    dlist = os.listdir( fdir )
    betti_list = [ f for f in dlist if f.endswith( suffix ) ]
    betti_list.sort(key=natural_key)

    print "chomp output list", betti_list
    
    # keep the frame numbers organized in a dict ?
    #betti = {}
    # nah, just list them
    betti_arr = []
    for b in betti_list:
        bnums = numpy.loadtxt( fdir+b, dtype=numpy.uint8 )
        betti_arr.append( bnums )
    betti_arr = numpy.asarray( betti_arr )
    return betti_arr.T
    
def plot_betti( barr, cell=1, savedir=None, dim=0, fig=None,
               total_cells=2, color='b' ):
    """
    Plot betti numbers for each frame for a cell. Obtain a time series
    (time=frame number)
    """
    if fig is None:
        fig = plt.figure()
    ax = fig.gca()
    #ax = fig.add_subplot(total_cells, 1, cell_num+1)
    data = barr[:,dim,:]
    ax.plot( data[1], 'o-', color=color, lw=1.5, ms=2 )
    # record title and some stats
    ax.set_title(  'Betti numbers for cell '+str(cell)+\
                 ' (mean='+str( round(data[1].mean()) )+')' )
    ax.set_xlabel( 'Frame' )
    ax.set_ylabel( r'$H_{'+str(dim)+'}$' )
    if savedir == None:
        fname = './figures_raw/betti_frames_H'+str(dim)+'_cell'+str(cell)+'.png'
    else:
        fname = savedir + '/betti_frames_H'+str(dim)+'_cell'+str(cell)+'.png'
    fig.savefig( fname )

def plot_hist( data, cell_num=1 ):

    fig = plt.figure()
    ax = fig.add_subplot(111)
    dmean = round( data.mean(), 1 )
    dvar = round( data.var(), 1 )
    n, bins, patches = ax.hist( data, 10, normed=1, facecolor='green',
                                alpha=0.75,  label=r"Mean="+str(dmean)+"\nVar="+str(dvar) )
    
    ## # add a 'best fit' line
    ## y = mlab.normpdf( bins, mu, sigma)
    ## l = plt.plot(bins, y, 'r--', linewidth=1)

    xmin = data.min()
    xmax = data.max()
    ax.set_xlabel(r'$H_1$ generators')
    ax.set_ylabel('Probability')
    ax.set_title(r"Distribution of $H_1$ Generators, Cell "+str(cell_num) )
    ax.axis([xmin-5, xmax+5, 0, 0.1])
    ax.grid(True)
    ax.legend()
    return fig

def plot_spectrum( data ):
    """
    Plot the power spectrum of the data.
    """
    d = data[1]
    # rfft gives positive frequecies. Square to get power spectrum.
    fp = numpy.absolute( numpy.fft.rfft( d ) )**2
    freq = numpy.fft.fftfreq( d.shape[-1] )
    n = len(fp)

    # reshape stuff a bit. keep only positive freqs.
    fp = fp[1:-1]
    freq = freq[1:n-1]
    lrslope = linregress( numpy.log(freq[30:]), numpy.log(fp[30:]) )[0]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.loglog( freq, fp, label="Lin. reg.="+str(round( lrslope,1 )) )
    ax.legend( loc='lower left' )
    return fig

###
# parallel functions
###
def chomp_stack( low, high, height, path, chomp_path, prefix ):

    for base in range( low,  high, height ):
        frames = []
        # list of frames to stack
        for x in range( height ):
            num = base + x
            frames.append( path + prefix + str( num ) + '.bmp' )
        stack = stack_images( frames, height )
        cubfile = chomp_path + prefix[:-1] + \
            '_b' + str( base ) + \
            '_h' + str( height ) 
        
        # Convert bmp files to array, stack them, write them to
        # chomp-readable format.
        array2chomp( stack, cubfile + '.cub' )
        
        # Now compute homology for each block
        run_chomp( cubfile + '.cub', cubfile + '.hom'  )

        
def run_mse( fname, **args ):
    """
    Input:
    -----

    fname : path to 1D data file containing time series to
    analyze. File can be .npy or .txt.

    If .npy, then file is loaded and rewritten to disk as a tmp file
    for mse input. Otherwise, if .txt (single column!), fname is
    passed straight to mse.

    **args:
    -----

    Eventually will correspond to MSE args (type 'mse --help' for more
    info).
    """
    try:
        data = numpy.load( fname )
        fname = fname[:-3] + 'txt'
        numpy.savetxt( fname, data )       
    except IOError:
        pass
        
    # form the command to pass to subprocess
    outfile = fname[:-3] + 'mse'
    cmd = 'mse <'+ fname + '>'+outfile
    try:
        p = subprocess.Popen( cmd, shell=True )
    except:
        print "Problems calling MSE"
        raise

def mse_converter( fname ):
    """
    Convert data in .mse file to a more friendly format. 
    """
    lines = []
    with open( fname ) as fh:
        for line in fh.readlines():
            if len(line) > 1:  # avoid empty lines
                if line.startswith('m'):
                    continue
                # strip off \n and split on tabs
                line = line.strip().split( '\t' )
                lines.append( ( float(line[0]), float(line[1]) ) )
    return numpy.array( lines )

        
if __name__ == "__main__":
    
    import pp
    import time

    start = time.time()

    stack_height = [10, 20, 30]

#    path = '/data/CT_Firn_Sample/output23-10-3/'
#    chomp_path = '/data/CT_Firn_Sample/chomp_files/'
#    prefix = 'K09b-23-10-'

    
    chomp_path = '/sciclone/data10/jberwald/CT_Firn_Samples/chomp_files/'
    path = '/sciclone/data10/jberwald/CT_Firn_Samples/output23-10-3/'
    prefix = 'K09b-23-10-'
    
    ## Write Chomp-readable files for 3D blocks

    #parallelize this stuff
#    ncpus = len( stack_height )
#    job_server = pp.Server( ncpus, ppservers=() )   
#    pool = []

#    bottom = 3200
#    top = 3220
    
    
    if 0:
        for height in stack_height:

            # pool.append( job_server.submit( chomp_stack,
            #                                 ( bottom,
            #                                   top,
            #                                   height,
            #                                   path,
            #                                   chomp_path,
            #                                   prefix ),
            #                                 depfuncs = ( stack_images, array2chomp,
            #                                              run_chomp )
            #                                 ) )

         
            print "Stack height:", height
            for base in range( 3500,  3700, height ):
                frames = []
                # list of frames to stack
                for x in range( height ):
                    num = base + x
                    frames.append( path + prefix + str( num ) + '.bmp' )

                print "    Stacking from base:", base
                stack = stack_images( frames, height )
                cubfile = chomp_path + prefix[:-1] + \
                    '_b' + str( base ) + \
                    '_h' + str( height ) 

                # Convert bmp files to array, stack them, write them to
                # chomp-readable format.
                array2chomp( stack, cubfile + '.cub' )

                # Now compute homology for each block
                run_chomp( cubfile + '.cub', cubfile + '.hom'  )
                extract_betti( cubfile + '.hom' )
            print ""

        print "Time:", time.time() - start

    if 1:
        for height in stack_height:
#        height = 10
            low = 3200
            high = 3400
#            low = 3500
#            high = 3700
            for dim in [0,1,2]:
                bettis = []
                for base in range( low,  high, height ):
                    betti_file = chomp_path + prefix[:-1] + \
                        '_b' + str( base ) + \
                        '_h' + str( height ) + \
                        '.betti'
                    bnums = numpy.loadtxt( betti_file, dtype=numpy.uint8 )
                    bettis.append( bnums[dim][1] )

                B = numpy.asarray( bettis, dtype=int )

                fig = plt.figure()
                ax = fig.gca()
                ax.plot( B, 'bo' )
                # draw a dashed line at the mean, same xlimits as below
                ax.hlines( B.mean(), -1, len(B), linestyle='dashed', colors='g', linewidth=2 )

                ax.set_title( r'$\beta_{'+str(dim)+'}$ for blocks of height='+\
                                  str(height)+'\nbetween frame '+str(low)+\
                                  ' and '+str(high) )
                ax.set_xlabel( "Block number (height="+str(height)+")" )
                ax.set_ylabel( r"$\beta_{"+str( dim )+"}$" )
                ax.set_xlim( -1, len(B) )
                ax.set_ylim( B.min() - 1, B.max() + 1 )

                fig_prefix = '/sciclone/data10/jberwald/CT_Firn_Samples/figures/'
                figname = 'binary_'+str( low )+ '_' + str( high ) +\
                    '_h'+str( height )+ '_d'+str( dim ) + '.png'
                fig.savefig( fig_prefix + figname, transparent=True )
                plt.close( fig )

