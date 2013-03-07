
import numpy
import os, re
import matplotlib.pyplot as plt
import cPickle as pkl
from pyTopTools.rbc_analysis import rbc_histogram as rh


def get_gens ( file, rmv='',data = ''):
    """
        Get generators function
    """
    with open(file, 'r') as f:
        s = f.read()
    f.close()
    #split up generators
    stringGens = s.split('\n')
    stringGens.remove('')
    gens = [map(int,sgen.split(' ')) for sgen in stringGens]
    if data:
        maxHeight = get_Max( data, 1)
        for gen in gens:
            if gen[-1] == -1:
                gen[-1] == maxHeight
    if rmv:
        remove_inf( gens )
    return gens

def get_gens_above( fname,epsilon ):
    """
        Epsilon1 is lower bound, epsilon2 is upper bound, i.e.
        epsilon1 < (death - birth) < epsilon2
    """
    gens = numpy.loadtxt( fname )
    goodGens = []
    for birth, death in gens:
        if ( death - birth ) > epsilon:
            goodGens.append( (birth,death) )
    return goodGens

def get_gens_between (file,epsilon1, epsilon2):
    """
        Epsilon1 is lower bound, epsilon2 is upper bound, i.e.
        epsilon1 < (death - birth) < epsilon2
    """
    with open(file, 'r') as f:
        s = f.read()
    f.close()
    goodGens = []
    #split up generators
    stringGens = s.split('\n')
    stringGens.remove('')
    gens = []
    #parse generators
    for sgen in stringGens:
        gens.append(map(int,sgen.split(' ')))
    for (birth,death) in gens:
        if (death - birth) > epsilon1 and (death-birth) < epsilon2:
            goodGens.append((birth,death))
    return goodGens

def get_gens_between_normed( fname, eps1, eps2, means=False ):
    """
    """
    gens = numpy.loadtxt( fname )
    y1, y2 = normalize_mid_lifespan( gens, eps1, eps2 )
        
    # normalize the generator stack and the midrange band
    gens = normalize( gens )
    goodGens = []
    # now find the normalized midrange generators birth and death time
    for (birth,death) in gens:
        if (death - birth) > y1 and (death-birth) < y2:
            goodGens.append((birth,death))
    if means:
        if goodGens:
            return numpy.asarray( goodGens )
        else:
            return None
    else:
        return goodGens

def normalize_mid_lifespan( gens, eps0, eps1 ):
    """
    Normalize the midrange band. Basically, take

    f: (eps0,eps1) --> [0,1].

    Returns f(eps0,eps1) = (y0,y1)
    """
    eps0 = float( eps0 )
    if eps1 == -1:
        eps1 = float( gens.max() )
    else:
        eps1 = float( eps1 )
    delta = gens.max()
    return eps0 / delta, eps1 / delta


def normalize(arr, imin=0, imax=1, dmin=None, dmax=None):
    """
    Normalize 'arr', in-place. (Stolen from stack
    overload. Surprised numpy doesn't have a built-in normalize
    function.)

    (imin, max) -- desired range of normalization

    dmin and dmax -- used if the array does not include all of the
    values. For example, birth time may not include the minimum and
    maximum values. In this case, 0 and max_height are passed to the
    function.
    """
    arr = arr.astype( float )
    if dmin is None:
        dmin = 0 #arr.min()
    if dmax is None:
        dmax = arr.max()
    arr -= dmin
    arr *= (imax - imin)
    arr /= (dmax - dmin)
    arr += imin
    return arr


def get_birth_times( cell, eps1, eps2=150, normed=True, first_bt=False ):
    """
    cell -- path to directory containing Perseus generator
    file (one for each frame).

    eps1 -- minimum lifespan

    eps2 -- maximum lifespan

    (So [eps1,eps2], stretched along the diagonal, is the band that we
    store generators from.)

    For each frame in cell==cell_dir: find generator lifespans.  Find first
    occurence, \tau, of a midrange generator ( use get_gens_between()
    for this, then peel off the birth time from the first
    (birth,death) pair to get birth time) store \tau.

    NOTE: Depending on eps1, very rarely a list of gens will be
    empty. We treat this as missing data and continue to loop over the
    frames.
    """
    frames = dir_list( cell )
    birth_times = []
    for frame in frames:
        if normed:
            gens = get_gens_between_normed( frame, eps1, eps2 )
        else:
            gens = get_gens_between( frame, eps1, eps2 )
        if not gens:
            continue
        if first_bt:
            birth_times.append( gens[0][0] )
        else:
            birth_times.append( gens )
    return birth_times


def get_midrange_gens( cell, eps1, eps2=-1, normed=False):
    """
    cell -- path to directory containing Perseus generator
    file (one for each frame).

    eps1 -- minimum lifespan

    eps2 -- maximum lifespan

    (So [eps1,eps2], stretched along the diagonal, is the band that we
    store generators from.)

    For each frame in cell==cell_dir: find generator lifespans.  Find first
    occurence, \tau, of a midrange generator ( use get_gens_between()
    for this, then peel off the birth time from the first
    (birth,death) pair to get birth time) store \tau.

    NOTE: Depending on eps1, very rarely a list of gens will be
    empty. We treat this as missing data and continue to loop over the
    frames.
    """
    frames = dir_list( cell )
    cell_stats = []
    for frame in frames[:10]:
        # get the midrange gen stats for frame (normed or not)
        if normed:
            gstats = get_gens_between_normed( frame, eps1, eps2, means=True )
        elif eps2 == -1:
            gstats = get_gens_above( frame, eps1 )
        else:
            gstats = get_gens_between( frame, eps1, eps2 )
        if not gstats:
            continue
        cell_stats.append( gstats )
    return cell_stats


def dir_list( fdir, betti=1 ):
    """
    Returns a list of Perseus output files for given betti #.
    """
    dlist = os.listdir( fdir )
    theFiles = [ fdir+f for f in dlist if f.endswith( '_'+str(betti)+'.txt' ) ]
    theFiles.sort( key=natural_key )
    return theFiles


def load_birth_times( old_hist = 'data/old_birth_times.pkl',
                      new_hist = 'data/new_birth_times.pkl' ):
    with open( old_hist ) as fh:
        old_bt = pkl.load( fh )
    with open( new_hist ) as fh:
        new_bt = pkl.load( fh )
    return old_bt, new_bt


def load_all_bt( prefix='' ):
    BT = defaultdict( dict )
    for eps in [30,40,50]:
        old_ = prefix+'./bt_old_normed_eps'+str(eps)+'.pkl'
        new_ = prefix+'./bt_new_normed_eps'+str(eps)+'.pkl'
        old, new = load_birth_times( old_, new_ )
        BT[eps]['old'] = old
        BT[eps]['new'] = new
    return BT


def natural_key(string_):
    """
    Use with frames.sort(key=natural_key)
    """
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string_)]
    
def characteristic_birth_time( fname, lower_bound, gen_num=1 ):
    """
    fname -- full path to file containing birth times for robust generators. Eg., old_gen_ts_lb40.pkl

    lower_bound -- the lower bound (eg., 40)
    
    gen_num -- 0 for first robust generator (probably the 'infinite' one), 1 for the first 'robust' generator, etc.

    """
    print "Getting birth times for ", fname 

    with open( fname ) as fh:
        cells = pkl.load( fh )
    cell_vals = cells.values()

    btimes = []
    # for each cell, stack the generators into one array
    for cell in cell_vals:
        for i, frame in enumerate( cell ):
            births = list( frame[:,0] )
            births.sort()
            try:
                btimes.append( births[gen_num] )
            except IndexError:
                print "bad birth list:", i, births
                continue
    return btimes 
                
def find_cell_avg( cell_concat, skip=3000 ):
    """
    Compute the average pixel height for a cell. 
    """
    cell = numpy.loadtxt( cell_concat, skiprows=skip )
    return cell.mean() 
    


if __name__ == "__main__":

    lb = [40,45,50,55,60,65,70,75]

    # CONSTRUCT ALL BIRTH TIMES OF GENERATORS ABOVE GIVEN THRESHOLD VALUE IN lb,
    # or compute averages
    if 1:
        new_prefix = '/data/PerseusData/PerseusOutput/original/2d_sparse/New/'
        old_prefix = '/data/PerseusData/PerseusOutput/original/2d_sparse/Old/'
        newlist = ['new_10', 'new_110125', 'new_130125', 'new_140125', 'new_3',
                   'new_4', 'new_40125', 'new_50125', 'new_6', 'new_60125', 'new_9']
        oldlist = ['old_100125', 'old_120125', 'old_15', 'old_2', 'old_4000', 'old_4001',
                   'old_5',  'old_50125',  'old_6',  'old_7',  'old_8',  'old_9',  'old_90125']

        new_avg_prefix = '/data/jberwald/wyss/data/Cells_Jesse/New/cells/'
        old_avg_prefix = '/data/jberwald/wyss/data/Cells_Jesse/Old/cells/'
        suffix = '-concatenated-ASCII'
        
        new_cells = [ prefix + c + '/' for c in newlist ]
        old_cells = [ prefix + c + '/' for c in oldlist ]

        old_ts = {}
        for x in lb:
            for cell in oldlist:
                print "Getting midrange gens for ", cell
                avg = find_cell_avg( old_avg_prefix + cell )
                bt = get_midrange_gens( old_prefix + cell+'/', x )
                old_ts[ cell ] = ( avg, bt )
            with open( '/data/jberwald/wyss/data/timeseries/rbc/old_gen_avg_lb'+str(x)+'.pkl', 'w' ) as fh:
                pkl.dump( old_ts, fh )

        new_ts = {}
        for x in lb:
            for cell in oldlist:
                print "Getting midrange gens for ", cell
                avg = find_cell_avg( new_avg_prefix + cell )
                bt = get_midrange_gens( new_prefix + cell +'/', x )
                new_ts[ cell ] = ( avg, bt )
            with open( '/data/jberwald/wyss/data/rbc/timeseries/new_gen_avg_lb'+str(x)+'.pkl', 'w' ) as fh:
                pkl.dump( new_ts, fh )

    if 0:
        
        prefix = '/sciclone/data10/jberwald/wyss/data/timeseries/'
        old_name = 'old_gen_ts_lb'
        new_name = 'new_gen_ts_lb'

        old_lb = {}
        for val in lb:
            with open( prefix + old_name + str( val ) + '.pkl' ) as fh:
                old = pkl.load( fh )
                v = old.values()
                # for each cell, stack the generators into one array
                allv = [ numpy.vstack( x ) for x in v ]
                # now stack all cell stacks into one array
                allv = numpy.vstack( allv )
            
            old_lb[ val ] = allv

        # dump all gens sorted by threshold to disk
        with open( prefix + 'old_orig_lb40-75.pkl', 'w' ) as fh:
            pkl.dump( old_lb, fh )

        new_lb = {}
        for val in lb:
            with open( prefix + new_name + str( val ) + '.pkl' ) as fh:
                new = pkl.load( fh )
                v = new.values()
                # for each cell, stack the generators into one array
                allv = [ numpy.vstack( x ) for x in v ]
                # now stack all cell stacks into one array
                allv = numpy.vstack( allv )
            
            new_lb[ val ] = allv
        
        # dump all gens sorted by threshold to disk
        with open( prefix + 'new_orig_lb40-75.pkl', 'w' ) as fh:
            pkl.dump( new_lb, fh )

    # EXTRACT THE K'TH BIRTH TIME
    if 0: 
        prefix = '/sciclone/data10/jberwald/wyss/data/timeseries/'
        old_name = 'old_gen_ts_lb'
        new_name = 'new_gen_ts_lb'

        gen_num = 1 # first birth time, or second,  or third, etc.

        for cell_type in [old_name, new_name]:
            for val in lb:
                bt = characteristic_birth_time( prefix + cell_type + str(val) + '.pkl',
                                                val,
                                                gen_num=gen_num )
                with open( prefix + cell_type + str(val) + '_gen'+str(gen_num)+'.pkl', 'w' ) as fh:
                    pkl.dump( bt, fh )
              

    if 0:
        new_cells = [ 'new_110125-concatenated-ASCII',
                      'new_140125-concatenated-ASCII',
                      'new_130125-concatenated-ASCII',
                      'new_40125-concatenated-ASCII',
                      'new_50125-concatenated-ASCII' ]
        old_cells = [ 'old_100125-concatenated-ASCII',
                      'old_50125-concatenated-ASCII',
                      'old_90125-concatenated-ASCII',
                      'old_120125-concatenated-ASCII']  
        prefix = '/sciclone/data10/jberwald/wyss/data/timeseries/'
        old_name = 'old_gen_ts_lb'
        new_name = 'new_gen_ts_lb'

        gen_num = 1 # first birth time, or second,  or third, etc.

        for cell_type in [old_name, new_name]:
            for val in lb:
                bt = characteristic_birth_time( prefix + cell_type + str(val) + '.pkl',
                                                val,
                                                gen_num=gen_num )
                with open( prefix + cell_type + str(val) + '_gen'+str(gen_num)+'.pkl', 'w' ) as fh:
                    pkl.dump( bt, fh )


    # MEANS
    if 0: 

        save_prefix = '/sciclone/data10/jberwald/wyss/data/timeseries/'
        prefix = '/sciclone/data10/jberwald/wyss/data/Cells_Jesse/'
        new_cells = [ 'new_110125-concatenated-ASCII',
                      'new_140125-concatenated-ASCII',
                      'new_130125-concatenated-ASCII',
                      'new_40125-concatenated-ASCII',
                      'new_50125-concatenated-ASCII' ]
        old_cells = [ 'old_100125-concatenated-ASCII',
                      'old_50125-concatenated-ASCII',
                      'old_90125-concatenated-ASCII',
                      'old_120125-concatenated-ASCII']
                      
        old_avgs = {}
        for cellname in old_cells:
            print "Getting average pixel height for ", cellname
            old_avgs[ cellname ] = find_cell_avg( prefix + 'Old/' + cellname )
        with open( save_prefix + 'old_avgs.pkl', 'w' ) as fh:
            pkl.dump( old_avgs, fh )

        new_avgs = {}
        for cellname in new_cells:
            print "Getting midrangeaverage pixel height for ", cellname
            new_avgs[ cellname ] = find_cell_avg( prefix + 'New/' + cellname )
        with open( save_prefix + 'new_avgs.pkl', 'w' ) as fh:
            pkl.dump( new_avgs, fh )

        print "AVERAGES:"
        print "old:", old_avgs
        print "new:", new_avgs

    # PLOT HISTOGRAMS
    if 0:
        #prefix = '/sciclone/data10/jberwald/wyss/data/timeseries/'
        prefix = '/data/jberwald/rbc/timeseries/'

        for val in lb:
            with open( prefix+ 'new_gen_ts_lb'+str(val)+'_gen1.pkl' ) as fh:
                new = pkl.load(fh)

            with open( prefix+ 'old_gen_ts_lb'+str(val)+'_gen1.pkl' ) as fh:
                old = pkl.load(fh)

            new = numpy.array( new, dtype=int )
            old = numpy.array( old, dtype=int )
            
            fig = rh.plot_hist( new, nbins=200 )
            fig = rh.plot_hist( old, nbins=200, fig=fig, color='r' )
            ax = fig.gca()
            ax.set_title( 'All birth time for lifespans above '+str(val) )
            plt.show()
            fig.savefig( './data/birth_times_lb'+str(val)+'_gen1.png' )
