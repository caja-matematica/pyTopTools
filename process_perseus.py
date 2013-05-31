import numpy as np
from pyTopTools.perseus_wrap import perseus
from pyTopTools.npy2perseus import write_cubtop, write_scubtop

def perseus_sub_super( mat, tmpPerseus, scale,
                       max_dim=1, pers_type='cubtop' ):
    """
    Performs matrix manipulation in order to compute persistent
    topology of superlevel sets as well as sublevel sets using Perseus
    software.
    
    population_mat -- An n x n matrix of values.

    tmpPerseus -- prefix of file name for writing to Perseus
    format. Should specify full path.

    scale -- Multiplicative factor by which to scale values in
    mat. Eg., since Perseus only deals with integers, this scale=100
    moves the decimal over 2 places. Hence, a greater value for scale
    give a greater resolution in the data.

    max_dim -- maximal dimension for combining Perseus output files.
    """

    # convert population matrix to integers
    M = scale * mat;
    M = np.ceil( M )
 
    # convert M to super-level format
    max_val = M.max() 
    Msup = -M + max_val + 1
    max_msup = Msup.max() 
    min_sup = Msup.min() 

    # Convert mat_text to Perseus file format. 
    perseus_in = tmpPerseus + '_sub.txt'
    perseus_in_sup = tmpPerseus + '_super.txt' 

    # assume if not sparse then full cubtop complex
    if pers_type == 'scubtop':
        print 'Writing sparse Perseus file!'
        write_scubtop( M, perseus_in );
        write_scubtop( Msup, perseus_in_sup )
    else:
        print 'Writing Perseus file!' 
        write_cubtop( M, perseus_in );
        write_cubtop( Msup, perseus_in_sup )
        
    # call assumes that 'perseus' is in your path. Perseus will
    # append '_*.txt' to perseus_out
    print 'Running Perseus on:' 
    print perseus_in 
    print perseus_in_sup 
    
    perseus_out = tmpPerseus + '_sub'
    perseus_out_sup = tmpPerseus + '_super'

    #run perseus on sub- and super-level versions of pop matrix
    perseus( perseus_in, perseus_out, pers_type )
    perseus( perseus_in_sup, perseus_out_sup, pers_type )
 
    # Now combine the output files into one file for diagram plotting
    for d in range( max_dim + 1 ):
     
        #print 'DIMENSION ' + str( d )
     
        tmp1 = perseus_out + '_' + str( d ) + '.txt'
        tmp2 = perseus_out_sup + '_' + str( d ) + '.txt'
        try:
            P = np.loadtxt( tmp1 );
            # if P is one line, then shape is (2,). reshape to (1,2)
            if len( P.shape ) == 1:
                P.resize( (1,2) )
                       
            # for dim > 0, sometimes we get [] matrices
            if P.size == 0:
                got_sub = False
            else:
                got_sub = True
        except IOError:
            got_sub = False
            print 'No perseus file found for ' + tmp1
        try:
            Psup = np.loadtxt( tmp2 );
            # if Psup is one line, then shape is (2,). reshape to
            # (1,2)
            if len( Psup.shape ) == 1:
                if len( Psup ) == 0:
                    got_super = False
                else:
                    Psup.resize( (1,2) )
            if Psup.size == 0:
                got_super = False
            else:
                got_super = True
        except IOError:
            got_super = False
            print 'No perseus file found for ' + tmp2
   
        # convert super-level Perseus data back using max_val and
        # account for -1's
        if got_super:
            # Shift births by max_val since we are considering
            # super-level sets.
            inf_idx = np.where( Psup[:,1] == -1 ) 
            non_inf = np.where( Psup[:,1] != -1 )[0]
            Psup = -Psup + max_val + 1
         
            # designate superlevel inf's differently and undo scaling
            Psup[ inf_idx, 1 ] = -2 
            Psup[ non_inf, 1 ] = Psup[ non_inf, 1 ] / scale
            Psup[:,0] = Psup[:,0] / scale
            if got_sub:
                non_inf = np.where( P[:,1] != -1 )[0]
                P[ non_inf, : ] = P[ non_inf, : ] / scale
                P[:,0] = P[:,0] / scale 
                full_pers = np.vstack( (P, Psup) )    
            else:
                # already rescaled above, so no need here
                full_pers = Psup

        # We just have sublevels 
        else:
            # don't need to deal with -1's
            non_inf = np.where( P[:,1] != -1 )[0] 
            P[ non_inf, 1 ] = P[ non_inf, 1 ] / scale
            P[:,0] = P[:,0] / scale 
            full_pers = P
     
        full_pers_fname = tmpPerseus + '_full_' + str( d ) + '.pdia'
        np.save( full_pers_fname, full_pers )
 
 
