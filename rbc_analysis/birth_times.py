
import numpy
import os
import matplotlib.pyplot as plt
import cPickle as pkl


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


if __name__ == "__main__":

    new_prefix = '/data/PerseusData/PerseusOutput/original/2d_sparse/New/'
    old_prefix = '/data/PerseusData/PerseusOutput/original/2d_sparse/Old/'
    newlist = ['new_10', 'new_110125', 'new_130125', 'new_140125', 'new_3',
               'new_4', 'new_40125', 'new_50125', 'new_6', 'new_60125', 'new_9']
    oldlist = ['old_100125', 'old_120125', 'old_15', 'old_2', 'old_4000', 'old_4001',
               'old_5',  'old_50125',  'old_6',  'old_7',  'old_8',  'old_9',  'old_90125']
    new_cells = [ prefix + c + '/' for c in newlist ]
    old_cells = [ prefix + c + '/' for c in oldlist ]
    frames = [ dir_list( c ) for c in cells ]
    # fig, ts = plot_hist_stack( frames, left_xlim=0.2, right_xlim=0.6, normed=False,
    #                            cutoff=0.2, ts_max=1000, skip=20, log=True )

    # lower bound
    lb = [40,45,50,55,60,65,70,75]
    old_ts = {}
    for x in lb:
        for cell_dir in old_cells:
            print "Getting midrange gens for ", cell_dir
            old_ts[ cell_dir ] = get_midrange_ts( cell_dir, x )
        with open( './timeseries/old_gen_ts_lb'+str(x)+'.pkl', 'w' ) as fh:
            pkl.dump( old_ts, fh )

    
    new_ts = {}
    for x in lb:
        for cell_dir in new_cells:
            print "Getting midrange gens for ", cell_dir
            new_ts[ cell_dir ] = get_midrange_ts( cell_dir, x )
        with open( './timeseries/new_gen_ts_lb'+str(x)+'.pkl', 'w' ) as fh:
            pkl.dump( new_ts, fh )
