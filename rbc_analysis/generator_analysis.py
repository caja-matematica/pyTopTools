import numpy as np 
import cPickle as pkl
from pyTopTools.rbc_analysis import birth_times as bt
from pyTopTools.rbc_analysis import robust_generators_hist as rg 
from matplotlib import pyplot as plt
import time

def load_gens( fname ):
    """
    fname : full path to file containing generators 
    """
    with open( fname, 'r' ) as fh:
        gens = pkl.load( fh )
    return gens

def threshold_frame( arr, eps ):
    """
    arr : N x 2 numpy array of (birth, death) pairs

    eps : return only death - birth > eps
    """
    idx = []
    for i, g in enumerate( arr ):
        if ( g[1] - g[0]) > eps:
            idx.append( i )
    return idx


def threshold_gens( gens, eps ):
    """gens : dictionary of generators, keyed by frame --> array of
    (birth,death) times.

    eps : lifespan cutoff. Eg., keep all (birth,death) pairs such that
    death-birth > eps.

    Returns dictionary keyed by frame number, which values containing
    arrays with only robust (birth,death) pairs.

    """
    robust = {}
    for frame_num, frame_gens in gens.iteritems():
        # indices of frame_gens with robust lifespans
        idx = threshold_frame( frame_gens, eps )
        robust[ frame_num ] = frame_gens[ idx ]
    return robust
    

def get_robust_generators( cell, eps=45, 
                           prefix='/sciclone/data10/jberwald/RBC/cells/persout/',
                           suffix='_gens.pkl' ):
    """cell : name of cell

    eps : lifespan threshold. Since we're just looking for the 'most
    robust' generators (top three, say), this is mostly a way to trim
    the noise off at this point.

    """
    cell_path = prefix + cell + suffix
    cell_gens = load_gens( cell_path )
    robust_gens = threshold_gens( cell_gens, eps )
    return robust_gens

def get_lifespans_all_frames( cell_gens ):
    """Concatenate lifespans over all frames for a single cell.

    cell : dictionary of cell birth,death generators keys by frame
    number.

    """
    lifespans = []
    for frame_num, frame_gens in cell_gens.iteritems():
        lifespans.append( frame_gens[:,1] - frame_gens[:,0] )
    # convert to list for consistency with original format
    return list( np.asarray( np.concatenate( lifespans ), dtype=np.int16 ) )

def save_gens( gens, fname ):
    """
    """
    with open( fname, 'w' ) as fh:
        pkl.dump( gens, fh )
        

def find_frame_avgs( cell_path, skip=0 ):
    """cell_path : full path to concatenated cell file

    skip : we don't really need all of the frames to get an idea of
    the average.

    """
    return bt.find_frame_avgs( cell_path, skip )

def find_cell_max( cell_path ):
    return bt.find_cell_max( cell_path )

def find_cell_avg( cell_path, skip=0 ):
    """cell_path : full path to concatenated cell file

    skip : we don't really need all of the frames to get an idea of
    the average.

    """
    return bt.find_cell_avg( cell_path, skip=2000 )

def concat_robust_generators( robust_all, max_gen=4 ):
    """robust_all : path to dictionary keyed by cells. Each cell maps to
    a dictionary
    
    max_gen : Search for up to the max_gen'th robust generator.

    """
    robust_gens = {}
    for g in range( 1, max_gen ):
        robust_gens[ g ] = bt.characteristic_birth_time( )

def get_frame_diff_vector( cell, boundary, skip=0 ):
    """
    Returns ( ||f_i - f_{i+1}||_{\infty} )
    """
    data = np.loadtxt( cell, skiprows=skip )
    bnd = np.loadtxt( boundary ).ravel()
    
    dvec = []
    # all operations are element-wise, so arrays are fine
    for i in xrange( len(data)-1 ):
        f0 = data[i] * bnd.ravel()
        f1 = data[i+1] * bnd.ravel()
        d = np.absolute( f0 - f1 )
        delta = d.max()
        dvec.append( delta )

    return dvec 
        
    
def load_pkls( new_name, old_name ):
    # load pickled data 
    with open( new_name ) as fh:
        new_data = pkl.load( fh ) 
    with open( old_name ) as fh:
        old_data = pkl.load( fh ) 
    return new_data, old_data


def plot_lags( new_fname, old_fname ):
    """Plot the norms of the lag vectors contained in the four file names
    above.

    On the x-axis is plotted the unnormalized data, on the y-axis the
    normalized data.

    """
    # extract the norms
    with open( new_fname + '.pkl' ) as fh:
        new_orig = pkl.load( fh )
    with open( new_fname + '_norm.pkl' ) as fh:
        new_norm = pkl.load( fh )
    new_orig = new_orig.values()
    new_norm = new_norm.values()

    # now the old cells
    with open( old_fname + '.pkl' ) as fh:
        old_orig = pkl.load( fh )
    with open( old_fname + '_norm.pkl' ) as fh:
        old_norm = pkl.load( fh )
    old_orig = old_orig.values()
    old_norm = old_norm.values()

    fig = plt.figure()
    ax = fig.gca()
    ax.plot( new_norm, new_orig, 'bo', ms=8 )
    ax.plot( old_norm, old_orig, 'ro', ms=8 )
    return fig

def plot_mean_number_of_generators( all_gens, eps, cell_type=None, fig=None, skip=1 ):
    """
    gens : array with (frame #, # robust gens in frame)

    eps : The epsilon threshold, used for plot legend/title

    fig : Optional figure object to pass in.   
    """
    if fig is None:
        fig = plt.figure( figsize=(10,8) )

    # find shortest block of frames
    min_gens = 5000
    for c,r in all_gens.items():
        if len( r ) < min_gens:
            min_gens = len( r )

    old = [ r[:min_gens,1] for c,r in all_gens.items() 
            if c.startswith('old') ]
    new = [ r[:min_gens,1] for c,r in all_gens.items() 
            if c.startswith('new') ]
    old = np.vstack( old )
    new = np.vstack( new )

    mean_old = old.mean( axis=0 )
    mean_new = new.mean( axis=0 )

    # Plot stuff
    ax = fig.gca()
    ax.plot( mean_old[::skip], 'r-', lw=2, 
             label='Old cells, $\epsilon=$'+str(eps) )
    
    ax.plot( mean_new[::skip], 'b-', lw=2, 
             label='Young cells, $\epsilon=$'+str(eps) )
    ax.set_xlabel( 'Frame number', fontsize=16 )
    ax.set_ylabel( 'Mean number of robust generators', fontsize=16 )
    ax.legend()
    fig.show()
    fig.savefig( 'figs_26cells/mean_robust_gens_eps'+str(eps)+'.pdf' )
    return fig
    

def get_number_of_robust_gens( gen_file ):
    """ number of robust generators. For each frame, subtract one to
    account for the 'infinite' generator.
    
    gen_file : pickled dictionary containing 

    { cells : { frame i : generator array } }
    """
    print "loading generator file..."
    with open( gen_file ) as fh:
        cells = pkl.load( fh )
    
    cell_info = {}
    for cell, frames in cells.iteritems():
        print cell
        gens_per_frame = []
        for i, gens in frames.iteritems():
            gens_per_frame.append( [ i, len( gens )-1 ] )
        cell_info[ cell ] = np.array( gens_per_frame )
    return cell_info
    



if __name__ == "__main__":

    
    prefix='/sciclone/data10/jberwald/RBC/cells/persout/'
    
    cells = []
    with open( 'cellnames' ) as fh:
        for line in fh.readlines():
            cells.append( line.strip() )

    new_cells = [ cell for cell in cells if cell.startswith('new') ]
    old_cells = [ cell for cell in cells if cell.startswith('old') ]
    # in-place reversal for switching cell pairings for comparison
    #old_cells.reverse()
    all_cells = zip( old_cells, new_cells )

    # compute the robust generators. The lifespan threshold eps==45
    # (default)
    if 0:
        robust_generators = {}
        for eps in [40,50,60,70]:
            print "eps", eps
            for cell in cells:
                print cell
                robust_generators[ cell ] = get_robust_generators( cell, eps=eps )
                
            save_gens( robust_generators, prefix+'robust_eps'+str(eps)+'.pkl' )  
            print ""

    # list of lists of each cell's lifespans (all). old and new separate
    if 0:
        new_ts = []
        old_ts = []
        for cell in cells:
            print "cell:", cell
            if cell.startswith( 'new' ):
                gens = load_gens( prefix + cell + '_gens.pkl' )
                # append the list of all lifespans over all frames of
                # a single cell
                new_ts.append( get_lifespans_all_frames( gens ) )
            else:
                gens = load_gens( prefix + cell + '_gens.pkl' )
                old_ts.append( get_lifespans_all_frames( gens ) )
        with open( prefix + 'new_hist_ts.pkl', 'w' ) as fh:
            pkl.dump( new_ts, fh )
        with open( prefix + 'old_hist_ts.pkl', 'w' ) as fh:
            pkl.dump( old_ts, fh )

    # statistics for cells and frames
    if 0:
        if 0:
            # first, find the avg pixel height for frame of each cell
            new_avgs = {}
            old_avgs = {}
            concat_prefix = '/sciclone/data10/jberwald/RBC/cells/'
            for cell in cells:
                print "cell:", cell
                cellname = cell + '-concatenated-ASCII'
                if cell.startswith( 'new' ):
                    # returns array of average for indexed by frame
                    new_avgs[ cell ] = find_frame_avgs( concat_prefix +'new/'+ cellname )
                else:
                    old_avgs[ cell ] = find_frame_avgs( concat_prefix +'old/'+ cellname )
            with open( prefix + 'new_avgs.pkl', 'w' ) as fh:
                pkl.dump( new_avgs, fh )            
            with open( prefix + 'old_avgs.pkl', 'w' ) as fh:
                pkl.dump( old_avgs, fh )      

        if 0:
            # first, find the avg pixel height over all frames of each cell
            new_avgs = {}
            old_avgs = {}
            concat_prefix = '/sciclone/data10/jberwald/RBC/cells/'
            for cell in cells:
                print "cell:", cell
                cellname = cell + '-concatenated-ASCII'
                if cell.startswith( 'new' ):
                    # returns array of average for indexed by frame
                    new_avgs[ cell ] = find_cell_avg( concat_prefix +'new/'+ cellname )
                else:
                    old_avgs[ cell ] = find_cell_avg( concat_prefix +'old/'+ cellname )
            with open( prefix + 'new_cell_avgs.pkl', 'w' ) as fh:
                pkl.dump( new_avgs, fh )            
            with open( prefix + 'old_cell_avgs.pkl', 'w' ) as fh:
                pkl.dump( old_avgs, fh )

        if 0:
            # first, find the max pixel height over all frames of each cell
            new_maxes = {}
            old_maxes = {}
            concat_prefix = '/sciclone/data10/jberwald/RBC/cells/'
            for cell in cells:
                tstart = time.time()
                print "cell:", cell
                cellname = cell + '-concatenated-ASCII'
                if cell.startswith( 'new' ):
                    # returns array of average for indexed by frame
                    new_maxes[ cell ] = find_cell_max( concat_prefix +'new/'+ cellname )
                    print "new maxes:", new_maxes
                else:
                    old_maxes[ cell ] = find_cell_max( concat_prefix +'old/'+ cellname )
                    print "old maxes:", old_maxes
                print "time:", time.time() - tstart
            with open( prefix + 'new_cell_maxes.pkl', 'w' ) as fh:
                pkl.dump( new_maxes, fh )            
            with open( prefix + 'old_cell_maxes.pkl', 'w' ) as fh:
                pkl.dump( old_maxes, fh )
      

        # we just separate the new and old robust generators here for
        # modularity
        if 0:
            robust_data = prefix + 'robust_all.pkl'
            with open( robust_data ) as fh:
                all_robust = pkl.load( fh )
            
            # split into new and old
            new_robust = {}
            old_robust = {}
            for cell in all_robust:
                print "cell:", cell
                if cell.startswith('new'):
                    new_robust[ cell ] = all_robust[ cell ]
                else:
                    old_robust[ cell ] = all_robust[ cell ]
            
            with open( prefix + 'new_robust.pkl', 'w' ) as fh:
                pkl.dump( new_robust, fh )
            with open( prefix + 'old_robust.pkl', 'w' ) as fh:
                pkl.dump( old_robust, fh )

        # separate the generators by their lifespans (or "robustness")
        if 0:
            for g in [1,2,3,4]:
                # new cells
                char_birth = bt.characteristic_birth_time( prefix+'new_robust.pkl', 
                                                           gen_num=g )
                with open( prefix + 'new_robust_gen'+str(g)+'.pkl', 'w' ) as fh:
                    pkl.dump( char_birth, fh )
                # ol cells
                char_birth = bt.characteristic_birth_time( prefix+'old_robust.pkl', 
                                                           gen_num=g )
                with open( prefix + 'old_robust_gen'+str(g)+'.pkl', 'w' ) as fh:
                    pkl.dump( char_birth, fh )

        # Separate generators, but subtract cell averages from the
        # robust gens
        if 0:

            death = False
            birth = True
            lspans = False
            keep_top = False
            norm = True

            #new_fname = 'new_birth_times_frameavg_shift_gen'
            #new_fname = 'new_birth_times_normed_gen'
            new_fname = 'new_birth_times_normed_low_gen'

            #old_fname = 'old_birth_times_frameavg_shift_gen'
            #old_fname = 'old_birth_times_normed_gen'
            old_fname = 'old_birth_times_normed_low_gen'
            
            # averages over all frames for each cell
            #new_avgs = 'new_cell_avgs.pkl'
            #old_avgs = 'old_cell_avgs.pkl'

            # avgs overfor each frame in each cell
            new_avgs = 'new_avgs.pkl' 
            old_avgs = 'old_avgs.pkl'

            for gen_num in [1,2]:
                print "gen:", gen_num

                # NEW CELLS
                new_char_data = bt.characteristic_generator_data_shift( prefix+'new_robust.pkl', 
                                                                        gen_num,
                                                                        shift_cell=False,
                                                                        shift_avg=prefix + new_avgs,
                                                                        lspans=lspans,
                                                                        keep_top=keep_top,
                                                                        death_times=death,
                                                                        birth_times=birth,
                                                                        norm_new=prefix+'new_cell_maxes.pkl')
                # if lspans:
                #     L = new_char_data[0]  
                #     #L = np.concatenate( L )

                # shift each frame by the overall cell avg
                # with open( prefix + 'new_birth_times_avgcellshift_gen'+\
                #            str(gen_num)+'.pkl', 'w' ) as fh:
                #     pkl.dump( new_char_data, fh )

                # shift each frame by frame avg
                if death:
                    with open( prefix + 'new_death_times_frameavg_shift_gen'+\
                               str(gen_num)+'.pkl', 'w' ) as fh:
                        pkl.dump( new_char_data, fh )
                # store all lifespans, trim the infinite one later
                elif lspans:
                    if keep_top:
                        with open( prefix + 'new_robust_lifespans_top.pkl', 'w' ) as fh:
                            pkl.dump( new_char_data, fh )
                    else:
                        with open( prefix + 'new_robust_lifespans_gen'+\
                                   str(gen_num) +'.pkl', 'w' ) as fh:
                            pkl.dump( new_char_data, fh )
                elif birth:
                    with open( prefix + new_fname + str(gen_num)+'.pkl', 'w' ) as fh:
                        pkl.dump( new_char_data, fh )

                # OLD CELLS
                old_char_data = bt.characteristic_generator_data_shift( prefix+'old_robust.pkl', 
                                                                        gen_num,
                                                                        shift_cell=False,
                                                                        shift_avg=prefix + old_avgs,
                                                                        lspans=lspans,
                                                                        keep_top=keep_top,
                                                                        death_times=death,
                                                                        birth_times=birth,
                                                                        norm_old=prefix+'old_cell_maxes.pkl')

                # if lspans:
                #     L = old_char_data[0] 
                #     #L = np.concatenate( L )

                # with open( prefix + 'old_birth_times_avgshift_gen'+\
                #            str(gen_num)+'.pkl', 'w' ) as fh:
                #     pkl.dump( old_char_data, fh )

                # shift each frame by frame avg
                # shift each frame by frame avg
                if death:
                    with open( prefix + 'old_death_times_frameavg_shift_gen'+\
                               str(gen_num)+'.pkl', 'w' ) as fh:
                        pkl.dump( old_char_data, fh )
                # store all lifespans, trim the infinite one later
                elif lspans:
                    if keep_top:
                        with open( prefix + 'old_robust_lifespans_top.pkl', 'w' ) as fh:
                            pkl.dump( old_char_data, fh )
                    else:
                        with open( prefix + 'old_robust_lifespans_gen'+\
                                   str(gen_num)+'.pkl', 'w' ) as fh:
                            pkl.dump( old_char_data, fh )
                elif birth:
                    with open( prefix + old_fname+ str(gen_num)+'.pkl', 'w' ) as fh:
                        pkl.dump( old_char_data, fh )
                
                print ""
   

        # Separate generators, normalizing them and sorting them by
        # birth threshold, not lifespan as above
        if 1:

            print "Finding lowest birth thresholds..."

            #new_fname = 'new_birth_times_frameavg_shift_gen'
            #new_fname = 'new_birth_times_normed_gen'
            new_fname = 'new_birth_times_normed_low_gen'

            #old_fname = 'old_birth_times_frameavg_shift_gen'
            #old_fname = 'old_birth_times_normed_gen'
            old_fname = 'old_birth_times_normed_low_gen'
            
            # averages over all frames for each cell
            #new_avgs = 'new_cell_avgs.pkl'
            #old_avgs = 'old_cell_avgs.pkl'

            # avgs overfor each frame in each cell
            new_avgs = 'new_avgs.pkl' 
            old_avgs = 'old_avgs.pkl'

            for gen_num in [0,1,2]:
                print "gen:", gen_num

                # NEW CELLS
                new_char_data = bt.characteristic_generator_data_birth_order( prefix+'new_robust.pkl', 
                                                                              gen_num,
                                                                              norm_new=prefix+'new_cell_maxes.pkl')

                with open( prefix + new_fname+str(gen_num)+'.pkl', 'w' ) as fh:
                    pkl.dump( new_char_data, fh )


                # OLD CELLS
                old_char_data = bt.characteristic_generator_data_birth_order( prefix+'old_robust.pkl', 
                                                                              gen_num,
                                                                              norm_old=prefix+'old_cell_maxes.pkl')

                with open( prefix + old_fname+str(gen_num)+'.pkl', 'w' ) as fh:
                    pkl.dump( old_char_data, fh )
                    
                print ""
   


        # plot birth time histograms for each generator levels this
        # portion can be copied to rbc_paper_figs.py for plotting
        # purposes.
        if 0:
            new_avgs = 'new_avgs.pkl' 
            old_avgs = 'old_avgs.pkl'
  
            # paths to files
            # new_birth_times = prefix + 'new_birth_times_avgcellshift_gen'
            # old_birth_times = prefix + 'old_birth_times_avgcellshift_gen'
            # new_birth_times = prefix + 'new_birth_times_frameavg_shift_gen'
            # old_birth_times = prefix + 'old_birth_times_frameavg_shift_gen'

            # new_birth_times = prefix + 'new_robust_gen'
            # old_birth_times = prefix + 'old_robust_gen'

            # new_birth_times = prefix + 'new_death_times_frameavg_shift_gen'
            # old_birth_times = prefix + 'old_death_times_frameavg_shift_gen'

            
            # new_birth_times = prefix + 'new_robust_lifespans_gen'
            # old_birth_times = prefix + 'old_robust_lifespans_gen'

            new_birth_times = prefix + 'new_birth_times_normed_gen'
            old_birth_times = prefix + 'old_birth_times_normed_gen'

            new_single = prefix + 'new_robust.pkl'
            old_single = prefix + 'old_robust.pkl'
            
            save_prefix = '/sciclone/data10/jberwald/RBC/cells/persout/figures/'
            
            death = False
            lspans = False
            normed = True

            if normed:
                norm_old = prefix+'old_cell_maxes.pkl'
                norm_new = prefix+'new_cell_maxes.pkl'
            else:
                norm_old = None
                norm_new = None

            for g in [1,2]:
                print "gen:", g
                # kind of a random pairing. We can alter this if there
                # are two especially "nice" cells
                for old,new in [all_cells[0]]:
                    print "cells:", old, new
                    if lspans:

                        new = old = None

                        old_fname = old_birth_times + str(g) + '.pkl'
                        new_fname = new_birth_times + str(g) + '.pkl'
                    else:

                        # remove this if we are overlaying single
                        # cells
                        #new = old = None

                        old_fname = old_birth_times+str(g)+'.pkl'
                        new_fname = new_birth_times+str(g)+'.pkl'

                    # plot the histograms
                    rg.birth_times_hist2( old_fname, 
                                          new_fname,
                                          save_prefix, g, 
                                          old, new, 
                                          old_single, new_single,
                                          prefix + old_avgs,
                                          prefix + new_avgs,
                                          death=death,
                                          lspans=lspans,
                                          normed=normed,
                                          norm_new=norm_new,
                                          norm_old=norm_old)
                print ""


    # Analyze lag vectors. computed separately in
    # compute_wasserstein.py shell scripts.
    if 0:
        import compute_wasserstein as C
        
        pnorm = [1,2,'inf']
        

        if 0:
            lags = [1,5,50]
            middle_term = '_pdist_lag'

            for lag in lags:
                print "lag:", lag
                for p in pnorm:
                    print "  pnorm:", p
                    newvals = {}
                    oldvals = {}
                    # find the norm of lag vector for each cell and
                    # each lag and p-norm
                    for cell in cells:
                        print "  cell:", cell
                        lagfile = prefix + cell + middle_term + str(lag) 
                        d = C.lag_vector( lagfile+'.npy', p=p )
                        if cell.startswith( 'new' ):
                            newvals[cell] = d
                        else:
                            oldvals[cell] = d

                    newsave = prefix + 'new_lag'+str(lag) + \
                              '_p' + str(p) +'.pkl'
                    oldsave = prefix + 'old_lag'+str(lag) + \
                              '_p' + str(p) +'.pkl'
                    with open( newsave, 'w' ) as fh:
                        pkl.dump( newvals, fh )
                    with open( oldsave, 'w' ) as fh:
                        pkl.dump( oldvals, fh )

        # now plot stuff
        if 1:
            lags = [1,5,50]
            for lag in lags:
                for p in pnorm:

                    # read in the files 
                    new_fname = prefix + 'new_lag'+str(lag) + \
                                '_p' + str(p) 
                    old_fname = prefix + 'old_lag'+str(lag) + \
                                '_p' + str(p) 

                    fig = plot_lags( new_fname, old_fname )
                    fig.savefig( 'new_old_lag'+str(lag)+'_p'+str(p)+'.png' )
                    
        # compute length of lag vectors for normalized frames
        if 0:
            middle_term = '_pdist_lag' 
            lags = [1,5,50]
            for lag in lags:
                print "lag:", lag
                for p in pnorm:
                    print "  pnorm:", p
                    newvals = {}
                    oldvals = {}
                    # find the norm of lag vector for each cell and
                    # each lag and p-norm
                    for cell in cells:
                        print "  cell:", cell
                        lagfile = prefix + cell + middle_term + str(lag) + '_norm'
                        d = C.lag_vector( lagfile+'.npy', p=p )
                        if cell.startswith( 'new' ):
                            newvals[cell] = d
                        else:
                            oldvals[cell] = d

                    newsave = prefix + 'new_lag'+str(lag) + \
                              '_p' + str(p) +'_norm.pkl'
                    oldsave = prefix + 'old_lag'+str(lag) + \
                              '_p' + str(p) +'_norm.pkl'
                    with open( newsave, 'w' ) as fh:
                        pkl.dump( newvals, fh )
                    with open( oldsave, 'w' ) as fh:
                        pkl.dump( oldvals, fh )

        
    # Normalized birth times for each cell individually
    if 0:
        norm_old = prefix+'old_cell_maxes.pkl'
        norm_new = prefix+'new_cell_maxes.pkl'

        new_robust = prefix + 'new_robust.pkl'
        old_robust = prefix + 'old_robust.pkl'
        
        old_bts = {}
        new_bts = {}
        for g in [1,2]:
            print "gen:", g
            for old_single, new_single in all_cells:
                # bt_new = bt.characteristic_generator_data_shift( new_robust, g, single_cell=new_single, 
                #                                                  birth_times=True, norm_new=norm_new )
                # bt_old = bt.characteristic_generator_data_shift( old_robust, g, single_cell=old_single, 
                #                                                  birth_times=True, norm_old=norm_old )


                # np.savetxt( './data/'+new_single+'_birthtimes_normed_g'+str(g)+'.txt', bt_new )
                # np.savetxt( './data/'+old_single+'_birthtimes_normed_g'+str(g)+'.txt', bt_old )


                bt_new = bt.characteristic_generator_data_birth_order( new_robust, g, 
                                                                       norm_new=norm_new,
                                                                       single_cell=new_single)
                bt_old = bt.characteristic_generator_data_birth_order( old_robust, g, 
                                                                       norm_old=norm_old,
                                                                       single_cell=old_single)


                np.savetxt( './data/'+new_single+'_birthtimes_normed_low_g'+str(g)+'.txt', bt_new )
                np.savetxt( './data/'+old_single+'_birthtimes_normed_low_g'+str(g)+'.txt', bt_old )


                # old_bts[ old_single ] = bt_old
                # new_bts[ new_single ] = bt_new

            
                print ""

    if 0:
        prefix = './data/robust_eps'
        skip = 10

        for eps in [40,50,60,70]:
            G = get_number_of_robust_gens( prefix + str( eps ) + '.pkl' )
            fig = plot_mean_number_of_generators( G, eps, skip=skip )
            
    if 0:
        from pyTopTools import chomp_tools as C
        import time

        for eps in [40,50,60,70]:
            print "eps:", eps

            # with open( prefix + 'robust_eps' + str(eps) + '.pkl' ) as fh:
            #     robust = pkl.load( fh )
            
            genfile = prefix + 'robust_eps' + str(eps) + '.pkl'
            G = get_number_of_robust_gens( genfile )
                
            for cell, vals in G.items():
                fname = prefix + cell + '_robust_eps' + str(eps)
                np.savetxt( fname+'.txt', vals[:,1] )

                C.run_mse( fname + '.txt' )
                time.sleep( 0.5 ) # give it time to write to disk?
                d = C.mse_converter( fname + '.mse' )
                
                fig_dir = '/sciclone/home04/jberwald/github/local/caja-matematica/pyTopTools/rbc_analysis/figs_26cells/'
                fig = plt.figure()
                ax = fig.gca()
                if cell.startswith( 'new' ):
                    color = 'b'
                else:
                    color = 'r'
                ax.plot( d[:,0], d[:,1], color+'-', lw=2 )
                ax.set_xlabel( 'Granularity', fontsize=16 )
                ax.set_ylabel( 'Sample Entropy', fontsize=16 )
                fig.savefig( fig_dir + cell + '_MSE_numgens_eps' + str(eps) +'.pdf' )
                
            print ""

    if 1:
        from pyTopTools import chomp_tools as C
        import time

        for cell in cells:
            print cell
            fname = prefix + cell + '_robust_longest_lifespan'
            C.run_mse( fname + '.txt' )
            time.sleep( 1 ) # give it time to write to disk?
            # returns 2-column time series (granularity x mse value)
            d = C.mse_converter( fname + '.mse' )

            # there seems to be some latency in disk reads that forces
            # the converter to find the file but return it
            # empty. since 'by hand' seems fine, can only assume that
            # if we keep trying we'll find a 'filled' file.
            while d.size == 0:
                print "trying to convert again..."
                d = C.mse_converter( fname + '.mse' )

            fig_dir = '/sciclone/home04/jberwald/github/local/caja-matematica/pyTopTools/rbc_analysis/figs_26cells/'
            fig = plt.figure()
            ax = fig.gca()
            if cell.startswith( 'new' ):
                color = 'b'
            else:
                color = 'r'
            ax.plot( d[:,0], d[:,1], color+'-', lw=2 )
            ax.set_xlabel( 'Granularity', fontsize=16 )
            ax.set_ylabel( 'Sample Entropy', fontsize=16 )
            ymin, ymax = ax.get_ylim()
            ax.set_ylim( (0, ymax+0.1) )
            fig.savefig( fig_dir + cell + '_MSE_longest_lifespan.pdf' )
                
            print ""
