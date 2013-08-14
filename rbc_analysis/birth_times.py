import numpy
import os, re
import matplotlib.pyplot as plt
import cPickle as pkl
import rbc_histogram as rh

from collections import defaultdict

def get_gens ( fname, rmv='',data = ''):
    """
        Get generators function
    """
    with open(fname, 'r') as f:
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

def continuous_features( data, gen_num=1, bound=40,
                         norm_old=None, norm_new=None,
                         birth=False, death=False,
                         longest_lifespan=False,
                         single_cell=None ):
    """Track threshold values of robust features in a continuous
    manner. Assumption:

    ||f_i - f_{i+1}|| < \eps,

    eg., the intensity function changes continuously. 

    1. Sort robustness by lifespan as in characteristic* below.
    
    2. Determine whether the latest feature's threshold is within 'bound' of
    the last birth time. (continuity assumption)

    3. If so, add it to the list of feature thresholds. 
    
    3'. If not, consider similarly robust lifespans (up or down of not
    searching for most robust) and search for a feature, say birth
    time, withing 'bound' of the last threshold.

    data : path to data file

    gen_num : the generator number, ordered by lifespan, with the
    infinite one pruned. first=1, second=2, etc. (not 0-based)

    """
    g = gen_num # account for 'infinite'
    
    try:
        print " Reading file..."
        # data is a file name
        with open( data ) as fh:
            cells = pkl.load( fh )
        cell_vals = cells.values()
    except TypeError:
        # in this case we assume that data contain a list of a single
        # cell's generators in each frame
        cell_vals = [ data ]

    # read in maxes if we want to normalize values
    if norm_old is not None:
        with open( norm_old ) as fh:
            maxes = pkl.load( fh )
    elif norm_new is not None:
        with open( norm_new ) as fh:
            maxes = pkl.load( fh )
    else:
        maxes = None

    # hold all the feature information (birth time, death time, or
    # longest lifespan)
    # if longest_lifespan:
    features = defaultdict( list )
    #        features = []

    bad_frames = 0

    # for each cell, stack the generators into one array
    for cell, cell_frames in cells.items():
        print "cell", cell
        # if we're just analyzing the birth times of a single cell,
        # move on if it's not the right cell. used for overlaying
        # single cells on birth times histogram.
        if single_cell is not None:
            if cell != single_cell: 
                continue 

        # list of robust threshold arrays in this frame
        cell_vals = cell_frames.values()

        # normalize generator valus by max pixel height if we have
        # maxes
        if maxes is not None:
            cell_maxes = maxes[ cell ]
            
        last_robust = None
        no_gen = False
        flipping = []

        # each frame contains (n x 2) arrays of (birth,death) thresholds
        for i, frame in enumerate( cell_vals ):

            lifespans = numpy.diff( frame )
            L = lifespans.ravel() 

            # find 'infinite' generator
            lmax = L.max()
            idx = numpy.where( L < lmax )[0]

            # idx of longest lifespan first
            sort_idx = L.argsort()[::-1]
            
            # g should be >= 1 (0 reserved for inf. gen)
            try:
                the_gen = sort_idx[ g ]
            # sometimes a frame is funky and there are no lifespans
            # greater than the cutoff except the infinite one
            except IndexError:
                continue

            if i == 0:
                if birth:
                    robust = frame[ the_gen, 0 ]  
                    life = L[ the_gen ]
                    features[cell].append( robust )
                    last_robust = robust
                    last_life = life
                elif death:
                    robust = frame[ the_gen, 1 ] 
                    life = L[ the_gen ]    
                elif longest_lifespan:
                    robust = L[ the_gen ]
                                   

            if i > 0:
                #robust, life = match_generators( frame, L, last_robust, last_life )
                w = numpy.where( L[idx] > 2*bound )[0]
                if len(w) == 0:
                    robust = frame[:,0].min()
                else:
                    robust = frame[w,0].min()
                if robust > 1000
                features[cell].append( robust )
                last_robust = robust
                last_life = life
    return features
            
    #         if birth or death:
    #             if last_robust is not None:
    #                 # find a feature that appears to be a closer match to
    #                 # our last robust threshold if 

    #                 match_generators( robust, life, last_robust, last_life )
    #                 if abs( robust - last_robust ) > bound:
    #                     flipping.append( i )
    #                     tries = g
    #                     while tries < len( L ) and abs( robust - last_robust ) > bound:
    #                         try:
    #                             the_gen = sort_idx[ g+1 ] 
    #                             robust = frame[ the_gen, axis ]
    #                             tries += 1
    #                         # no generators fit the continuity criterion
    #                         except IndexError:
    #                             no_gen = True
    #                             break
            
    #         if no_gen:
    #             # print "NO GEN"
    #             # print "frame num:", i
    #             # print "  last:", last_robust, last_life
    #             # print "  current:", frame
    #             # print "  lifespans:", L
    #             # bad_frames += 1
    #             continue
    #         else:
    #             # normalize generator valus by max pixel for this cell
    #             if maxes is not None:
    #                 robust = float( robust ) / cell_maxes
    #             if birth or death:
    #                 features.append( robust )
    #             else:
    #                 features[ cell ].append( robust )
    #             last_robust = robust
    #             last_life = life
           
            
    #     print ""
    #     #print "BAD FRAMES", bad_frames
    # return features


def match_generators( current_frame, current_life, last_robust, last_life ):
    """Find the birth time in g2 that most closely matches that in g1 in
    terms of 1) threshold level and 2) lifespan.

    gen : first generator to try
    """
    life_diff = numpy.abs( current_life - last_life )
    gen_diff = numpy.abs( current_frame[:,0] - last_robust )

    life_idx = life_diff.argsort()
    gen_idx = gen_diff.argsort()

    # find the closest birth threshold
    if gen_idx[0] == life_idx[0]:
        closest_gen = current_frame[ gen_idx[0],0 ]
        life = current_life[ life_idx[0] ]
    elif gen_idx[1] == life_idx[1]:
        closest_gen = current_frame[ gen_idx[1],0 ]
        life = current_life[ life_idx[1] ]
    # elif gen_idx[0] == life_idx[1]:
    #     closest_gen = current_frame[ gen_idx[0],0 ]
    #     life = current_life[ life_idx[1] ]
    # elif gen_idx[1] == life_idx[0]:
    #     closest_gen = current_frame[ gen_idx[1],0 ]
    #     life = current_life[ life_idx[0] ]
    # just resort to grabbing the closest lifespan generator
    else:
        closest_gen = current_frame[ life_idx[0],0 ]
        life = current_life[ life_idx[0] ]
    
    return closest_gen, life

    # # where do things overlap
    # optimal = (life_idx == gen_idx)
    # matches = numpy.where( optimal )[0]

    # # best match
    # if len( matches ) > 0:
    #     new_robust = current_frame[ matches[0], 0 ]
    #     new_life = current_life[ matches[0] ]
    #     return new_robust, new_life
    # else:
    #     return last_robust, last_life

def characteristic_birth_time( data, gen_num=1, top=False, lspans=False ):
    """
    data -- full path to file containing birth times for robust
    generators or list of cell's generators. Eg., old_gen_ts_lb40.pkl

    gen_num -- 0 for first robust generator (probably the 'infinite'
    one), 1 for the first 'robust' generator, etc.

    top -- bool. If True, grab the top gen_num robust generators, as
    opposed to the default, which is to grab only the gen_num'th
    generator.
    
    lspans : bool. return the lifespan array as well

    """
    g = gen_num + 1 # account for 'infinite'
    
    try:
        # data is a file name
        with open( data ) as fh:
            cells = pkl.load( fh )
        cell_vals = cells.values()
    except TypeError:
        # in this case we assume that data contain a list of a single
        # cell's generators in each frame
        cell_vals = [ data ]

    btimes = []
    
    # save lifespans for testing separation of generators
    if lspans:
        LS = []
        
    # for each cell, stack the generators into one array
    for cell in cell_vals:
        # small hack for new list of cells (7/1/13)
        if type( cell ) == dict:                       
            cell = cell.values()  
        for i, frame in enumerate( cell ):
            lifespans = numpy.diff( frame )
            L = lifespans.T[0]
            if lspans:
                L.sort() # in-place sort
                LS.append( L )
            # argsort give indices of sorted array, from lowest to
            # highest (take only largest 'gen_num' lifespans; final
            # [::-1] reverses order). so the 'infinite' generator
            # lifespan is last. using -g gets us (gen_num+1) slots
            # back. When all is done, the birth time of the infinite
            # generator is at index 0.
            sort_idx = L.argsort()[-g:][::-1]
            births = list( frame[sort_idx,0] )
            births.sort() 
            try:
                if top:
                    # excludes infinite generator which will come
                    # first in the list.
                    btimes.append( tuple(births[1:g]) )
                else:
                    btimes.append( births[gen_num] )
            # sometimes there aren't enough robust generators found in
            # a frame for the gen_num index. The number per cell is
            # typically fairly small so we treat these as outliers.
            except IndexError:
                #print "bad birth list:", i, births
                continue
    if lspans:
        return btimes, LS
    else:
        return btimes

def characteristic_generator_data_birth_order( data, gen_num, 
                                               norm_old=None, norm_new=None,
                                               single_cell=None):
    """Shift the birth times by average contained in shift_avg. These are
    either per-cell or per-frame per-cell.


    data -- full path to file containing birth and death times for robust
    generators or list of cell's generators. Eg., old_gen_ts_lb40.pkl

    gen_num -- 0 for first robust generator (probably the 'infinite'
    one), 1 for the first 'robust' generator, etc.

    shift_avg : path to file containing cell avg or frame averages

    shift_cell : bool. compute the shift based on entire cell average
    or frame by frame.

    top -- bool. If True, grab the top gen_num robust generators, as
    opposed to the default, which is to grab only the gen_num'th
    generator.
    
    norm_{old,new} : path to frame max values for normalizing

    lspans : bool. return the lifespan array as well

    """
    g = gen_num # account for 'infinite'
    
    try:
        # data is a file name
        with open( data ) as fh:
            cells = pkl.load( fh )
        cell_vals = cells.values()
    except TypeError:
        # in this case we assume that data contain a list of a single
        # cell's generators in each frame
        cell_vals = [ data ]


    # read in maxes if we want to normalize values
    if norm_old is not None:
        with open( norm_old ) as fh:
            maxes = pkl.load( fh )
    elif norm_new is not None:
        with open( norm_new ) as fh:
            maxes = pkl.load( fh )
    else:
        maxes = None

    # hold all the birth times 
    btimes = []
          
    # for each cell, stack the generators into one array
    for cell, cell_frames in cells.items():
        print "cell", cell
        # if we're just analyzing the birth times of a single cell,
        # move on if it's not the right cell. used for overlaying
        # single cells on birth times histogram.
        if single_cell is not None:
            if cell != single_cell: 
                continue 

        # list of robust threshold arrays in this frame
        cell_vals = cell_frames.values()

        # normalize generator valus by max pixel height if we have
        # maxes
        if maxes is not None:
            cell_maxes = maxes[ cell ]
            
        # each frame contains (n x 2) arrays of (birth,death) thresholds
        for i, frame in enumerate( cell_vals ):

            # normalize generator valus by max pixel for this cell
            if maxes is not None:
                frame = numpy.asarray( frame, dtype=numpy.float )
                frame /= cell_maxes

            lifespans = numpy.diff( frame )
            L = lifespans.ravel() 

            # find 'infinite' generator
            lmax = L.max()
            idx = numpy.where( L < lmax )[0]
            
            # new array that excludes max
            robust = frame[idx,0]
            
            # in-place sort: n'th born robust generator in position
            # zero.
            robust.sort()

            # g = 0 (first born robust), g=1 (second born), etc.
            try:
                btimes.append( robust[g] )
            except IndexError:
                print "lifespan", L
                print "robust", robust
                print "g:", g
    
    return btimes

def characteristic_generator_data_shift( data, gen_num, 
                                         shift_cell=False,  shift_avg=None,
                                         shift_frames=False,
                                         keep_top=False, lspans=False,
                                         single_cell=True, death_times=False,
                                         birth_times=False,
                                         norm_old=None, norm_new=None ):
    """Shift the birth times by average contained in shift_avg. These are
    either per-cell or per-frame per-cell.


    data -- full path to file containing birth and death times for robust
    generators or list of cell's generators. Eg., old_gen_ts_lb40.pkl

    gen_num -- 0 for first robust generator (probably the 'infinite'
    one), 1 for the first 'robust' generator, etc.

    shift_avg : path to file containing cell avg or frame averages

    shift_cell : bool. compute the shift based on entire cell average
    or frame by frame.

    top -- bool. If True, grab the top gen_num robust generators, as
    opposed to the default, which is to grab only the gen_num'th
    generator.
    
    norm_{old,new} : path to frame max values for normalizing

    lspans : bool. return the lifespan array as well

    """
    g = gen_num + 1 # account for 'infinite'
    
    try:
        # data is a file name
        with open( data ) as fh:
            cells = pkl.load( fh )
        cell_vals = cells.values()
    except TypeError:
        # in this case we assume that data contain a list of a single
        # cell's generators in each frame
        cell_vals = [ data ]

    # each cell's average. must also provide shift_avg path
    if shift_cell:
        if shift_avg is None:
            print "Must provide path to averages"
        with open( shift_avg ) as fh:
            shift = pkl.load( fh )

    # read in maxes if we want to normalize values
    if norm_old is not None:
        with open( norm_old ) as fh:
            maxes = pkl.load( fh )
    elif norm_new is not None:
        with open( norm_new ) as fh:
            maxes = pkl.load( fh )
    else:
        maxes = None

    # hold all the birth times 
    btimes = []
    
    # save lifespans for testing separation of generators
    if lspans:
        LS = []
        
    # for each cell, stack the generators into one array
    for cell, cell_frames in cells.items():
        print "cell", cell
        # if we're just analyzing the birth times of a single cell,
        # move on if it's not the right cell. used for overlaying
        # single cells on birth times histogram.
        if single_cell is not None:
            if cell != single_cell: 
                continue 

        # list of robust threshold arrays in this frame
        cell_vals = cell_frames.values()

        # normalize generator valus by max pixel height if we have
        # maxes
        if maxes is not None:
            cell_maxes = maxes[ cell ]
            
        # each frame contains (n x 2) arrays of (birth,death) thresholds
        for i, frame in enumerate( cell_vals ):

            # normalize generator valus by max pixel for this cell
            if maxes is not None:
                frame = numpy.asarray( frame, dtype=numpy.float )
                frame /= cell_maxes

            lifespans = numpy.diff( frame )
            L = lifespans.ravel() #T[0]
            if lspans:
                L.sort() # in-place sort
                LS.append( L )

            # argsort give indices of sorted array, from lowest to
            # highest (take only largest 'gen_num' lifespans; final
            # [::-1] reverses order). so the 'infinite' generator
            # lifespan is last. using -g gets us (gen_num+1) slots
            # back. When all is done, the birth time of the infinite
            # generator is at index 0.
            sort_idx = L.argsort()[-g:][::-1]
    
            # assumes that each 'shift' is a dict
            # {cell:avg_val_over_all_frames}
            if shift_cell:
                if death_times:
                    coord = frame[sort_idx,1] - shift[cell]
                elif birth_times:
                    coord = frame[sort_idx,0] - shift[cell] 
                else:
                    print "No rule for how to shift cell"
                #coord = numpy.asarray( coord, dtype=numpy.int )
            # we'll have to shift frame by frame below
            elif shift_frames:
                if death_times:
                    coord = frame[sort_idx,1] - shift[cell][i]
                elif birth_times:
                    coord = frame[sort_idx,0] - shift[cell][i]
                else:
                    print "No rule for how to shift cell"
            # no shifting, just the values, please
            else:
                if death_times:
                    coord = frame[sort_idx,1]
                elif birth_times:
                    coord = frame[sort_idx,0]
                # reorder lifespans from largest to smallest
                elif lspans:
                    coord = L[sort_idx]
                else:
                    coord = frame[sort_idx]
                
            # if normed, casting as int will zero out values
            if maxes is None:
                coord = numpy.asarray( coord, dtype=numpy.int )
            # else:
            #     coord = numpy.asarray( coord, dtype=numpy.float )
                                
            coord = list( coord )
            #coord.sort() 
            try:
                if keep_top:
                    # excludes infinite generator which will come
                    # first in the list, but keeps all the gens, sorted.
                    btimes.append( coord[1:g] )
                else:
                    btimes.append( coord[gen_num] )
            # sometimes there aren't enough robust generators found in
            # a frame for the gen_num index. The number per cell is
            # typically fairly small so we treat these as outliers.
            except IndexError:
                #print "bad birth list:", i, coord
                continue
    # if lspans:
    #     return btimes, LS
    #else:
    return btimes
                
def find_cell_avg( cell_concat, skip=3000 ):
    """
    Compute the average pixel height for a cell. 
    """
    cell = numpy.loadtxt( cell_concat, skiprows=skip )
    #w = nump.where( cell!= 0 )
    # find the mean of all non-zero pixels. 6/29/13
    return cell.mean()

def find_frame_avgs( cell_concat, skip=0 ):
    """
    Find the average pixel value for each frame
    """
    cell = numpy.loadtxt( cell_concat, skiprows=skip )
    cell_avgs = []
    # grab frame number (row index) and frame (row)
    for frame in cell:
        w = numpy.where( frame != 0 )
        cell_avgs.append( frame[w].mean() )
    return numpy.asarray( cell_avgs )

def find_cell_max( cell_concat ):
    cell = numpy.loadtxt( cell_concat )
    return cell.max()

def find_frame_maxes( cell_concat ):
    cell = numpy.loadtxt( cell_concat )
    cell_maxes = []
    for frame in cell:
        w = numpy.where( frame != 0 )
        cell_maxes.append( frame[w].mean() )
    return numpy.asarray( cell_maxes )
    

def align_avg_birth_times( avgs, bts ):

    alignment = {}
    
    avg_names = []
    for x in avgs:
        xs = x.split('-')[0]
        avg_names.append( xs )

    for x in bts:
        xs = x.split('/')[-2]
        if xs in avg_names:
            alignment[ xs ] = bts[ x ]

    return alignment
    

def shift_birth_times( birth_times, avgs ):
    """
    birth_times -- birth times of robust generators for all frames for
    a single cell.

    shift -- the average pixel value for the cell giving birth_times,
    shift all generators by this amount as a way to normalize across cells.
    """
    pass

def save_longest_lifespans():
    """Convenience function to save text files of the 'longest lifespan'
    arrays for the each cell as a text file.
    """
    new_robust = '/sciclone/data10/jberwald/RBC/cells/persout/new_robust.pkl'
    old_robust = '/sciclone/data10/jberwald/RBC/cells/persout/old_robust.pkl'

    if norm:
        norm_old = prefix+'old_cell_maxes.pkl'
        norm_new = prefix+'new_cell_maxes.pkl'
    else:
        norm_old = None
        norm_new = None

    new = continuous_features( new_robust, 1, longest_lifespan=True, norm_new=norm_new )
    old = continuous_features( old_robust, 1, longest_lifespan=True, norm_old=norm_old )
    
    prefix = '/sciclone/data10/jberwald/RBC/cells/persout/'
    for cell, vals in new.items():
        numpy.savetxt( prefix + cell + '_robust_longest_lifespan.txt', vals )
        
    for cell, vals in old.items():
        numpy.savetxt( prefix + cell + '_robust_longest_lifespan.txt', vals )

def plot_longest_lifespans():
    """
    """
    new_robust = '/sciclone/data10/jberwald/RBC/cells/persout/new_robust.pkl'
    old_robust = '/sciclone/data10/jberwald/RBC/cells/persout/old_robust.pkl'
    new = continuous_features( new_robust, 1, longest_lifespan=True )
    old = continuous_features( old_robust, 1, longest_lifespan=True )

    save_prefix = '/sciclone/home04/jberwald/github/local/caja-matematica/pyTopTools/rbc_analysis/figs_26cells/'

    for cell, vals in new.items():
        print "Plotting ", cell 
        fig = plt.figure()
        ax = fig.gca()
        ax.plot( vals, 'b-', label=cell )
        ax.set_xlabel( 'time', fontsize=16 )
        ax.set_ylabel( 'Lifespan of longest-lived robust generator', fontsize=16 )
        ax.legend()
        fig.savefig( save_prefix + cell + '_longest_lifespan.pdf' )
        

    for cell, vals in old.items():
        print "Plotting ", cell
        fig = plt.figure()
        ax = fig.gca()
        ax.plot( vals, 'r-', label=cell )
        ax.set_xlabel( 'time', fontsize=16 )
        ax.set_ylabel( 'Lifespan of longest-lived robust generator', fontsize=16 )
        ax.legend()
        fig.savefig( save_prefix + cell + '_longest_lifespan.pdf' )
    



if __name__ == "__main__":

    lb = [40,45,50,55,60,65,70,75]

    # CONSTRUCT ALL BIRTH TIMES OF GENERATORS ABOVE GIVEN THRESHOLD VALUE IN lb,
    # or compute averages
    if 0:
        new_prefix = '/data/PerseusData/PerseusOutput/original/2d_sparse/New/'
        old_prefix = '/data/PerseusData/PerseusOutput/original/2d_sparse/Old/'
        newlist = ['new_10', 'new_110125', 'new_130125', 'new_140125', 'new_3',
                   'new_4', 'new_40125', 'new_50125', 'new_6', 'new_60125', 'new_9']
        oldlist = ['old_100125', 'old_120125', 'old_15', 'old_2', 'old_4000', 'old_4001',
                   'old_5', 'old_50125', 'old_6', 'old_7', 'old_8', 'old_9', 'old_90125']

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
            for cell in newlist:
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
              
    # GRAB GENERATOR NUMBER X 
    if 0:
        prefix = '/sciclone/data10/jberwald/wyss/data/timeseries/'
        old_name = 'oldgen_avgshift_lb'
        new_name = 'newgen_avgshift_lb'

        gen_nums = [1,2,3,4,5]  # first birth time, or second,  or third, etc.

        for cell_type in [old_name, new_name]:
            for val in lb:
                bt = characteristic_birth_time( prefix + cell_type + str(val) + '.pkl',
                                                val,
                                                gen_num=gen_num )
                with open( prefix + cell_type + str(val) + '_gen'+str(gen_num)+'.pkl', 'w' ) as fh:
                    pkl.dump( bt, fh )

        for g in gen_nums:
            for cell_type in [old_name, new_name]:
                for val in lb:
                    bt = characteristic_birth_time( prefix + cell_type + str(val) + '.pkl',
                                                    val, gen_num=g )
                    with open( prefix + cell_type + str(val) + '_gen'+str(g)+'.pkl', 'w' ) as fh:
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
    if 1:
        #prefix = '/sciclone/data10/jberwald/wyss/data/timeseries/'
        prefix = '/data/jberwald/rbc/timeseries/'

        oldname = 'oldgen_avgshift_lb'
        newname = 'newgen_avgshift_lb'

        new_single = 'newgen_avgshift_lb40.pkl'
        old_single = 'oldgen_avgshift_lb40.pkl'

        new_cells = [ 'new_50125', 'new_40125', 'new_130125' ]
        old_cells = [ 'old_100125', 'old_50125', 'old_90125' ]

        gens = [1,2]#,3,4]
        val = 40  # lower bound shouldn't matter, so just grab the lb=40 file

        for c_new, c_old in zip( new_cells, old_cells ):
            for g in gens:
                #for val in lb:
                with open( prefix+ newname +str(val)+'_gen'+str(g)+'_V2.pkl' ) as fh:
                    new = pkl.load(fh)

                with open( prefix+ oldname +str(val)+'_gen'+str(g)+'_V2.pkl' ) as fh:
                    old = pkl.load(fh)

                ####################
                ## FOR OVERLAYING SINGLE CELLS
                with open( prefix+ new_single ) as fh:
                    a = pkl.load(fh)
                    sing_new = a[ c_new ]

                with open( prefix+ old_single ) as fh:
                    b = pkl.load(fh)
                    sing_old = b[ c_old ]

                bt_new = characteristic_birth_time( sing_new, gen_num=g )
                bt_new = numpy.array( bt_new )
                bt_old = characteristic_birth_time( sing_old, gen_num=g )
                bt_old = numpy.array( bt_old )
                #
                ####################

                new = numpy.array( new, dtype=int )
                old = numpy.array( old, dtype=int )

                fig = rh.plot_hist( new, nbins=100 )
                fig = rh.plot_hist( old, nbins=100, fig=fig, color='r' )

                fig = rh.plot_hist( bt_new, nbins=100, fig=fig, color='g' )
                fig = rh.plot_hist( bt_old, nbins=100, fig=fig, color='y' )

                ax = fig.gca()
                # ax.set_title( r'Birth times for robust generator #'+str(g)+\
                #                 ' (lifespan threshold = '+str(val)+')', fontsize=16 )
                ax.set_xlabel( r'Birth threshold (' +c_new+' and '+ c_old +')', fontsize=16 )
                ax.set_ylabel( r'Number of generators ($g='+str(g)+'$)', fontsize=16 )
                plt.show()
                fig.savefig( './data/birth_times_gen'+str(g)+'_' + c_old +'_'+c_new + '.png' )
