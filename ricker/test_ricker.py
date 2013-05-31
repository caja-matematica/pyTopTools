#!/usr/bin/python
"""
Test persistence diagram structure against spatial correlation measures.

Tests:

checkerboard
    - +/-1 (scaled) and random values in [1,2] (scaled)

patches/islands
    - small --> large (same topology but different spatial correlation)

rings
    - small --> large (ditto)
"""
from pyTopTools.ricker import ricker_patch_model as rpm
import numpy as np

def test_checkerboard_binary( n, trials=1 ):
    """
    For a -1/+1 checkerboard, Moran's I should be -1. 
    """
    mat = 2 * np.ones( (n,n), dtype=int )
    mat = mat.ravel()
    nx = len( mat )

    mat.resize( (n,n) )
    for i in range( n ):
        #  if i % 2 == 0:
        print i
        mat[i,i] = 1
    mat[0,-1] = 1
    mat[-1,0] = 1

    #mat.resize( (n,n) )

    print mat
            
    C = rpm.Ricker( mat )

    # morans I
    mi = C.compute_moran()
    
    # PD
    # shift mean >= 1
    C.compute_persistence( 'debug/spatial_checkerBINARY_', 1, shift_mean=2 )
    fig0 = C.draw_persistence_diagram( 'debug/spatial_checkerBINARY_' +\
                                       '_full_'+str(0)+'.pdia.npy' )
    fig1 = C.draw_persistence_diagram( 'debug/spatial_checkerBINARY_' +\
                                       '_full_'+str(1)+'.pdia.npy' )
    fig0.savefig( 'debug/spatial_checkerBINARY_' + 'full_'+str(0)+'.png' )
    fig1.savefig( 'debug/spatial_checkerBINARY_' + 'full_'+str(1)+'.png' )
    figs = ( fig0, fig1 )
    return mi, figs, mat
 
def test_checkerboard( n, trials=1 ):
    """
    For a -1/+1 checkerboard, Moran's I should be -1. 
    """
    mat = np.empty( (n,n) )
    mat = mat.ravel()
    nx = len( mat )

    print mat

    morans = []
    pdfigs = []
    patches = []
    
    for trial in range( trials ):
        A = mat.copy()
        for i in range( nx ):
            if i % 2 == 0:
                # above zero
                A[i] = 2 + np.random.uniform()
            else:
                # below zero patches
                A[i] = 1 + np.random.uniform()
        A.resize( (n,n) )
        patches.append( A )
        C = rpm.Ricker( A )

        # morans I
        morans.append( C.compute_moran() )

        # PD
        C.compute_persistence( 'debug/spatial_checkerboard_'+str( trial ), scale=1000, shift_mean=2 )
        fig0 = C.draw_persistence_diagram( 'debug/spatial_checkerboard_'+str( trial ) +\
                                           '_full_'+str(0)+'.pdia.npy' ) 
        fig1 = C.draw_persistence_diagram( 'debug/spatial_checkerboard_'+str( trial ) +\
                                           '_full_'+str(1)+'.pdia.npy' )
        fig0.savefig( 'debug/spatial_checkerboard_trial'+str( trial ) + '_full_'+str(0)+'.png' )
        fig1.savefig( 'debug/spatial_checkerboard_trial'+str( trial ) + '_full_'+str(1)+'.png' )
        
        d0 = np.load( 'debug/spatial_checkerboard_'+str( trial ) +\
                      '_full_'+str(0)+'.pdia.npy' )
        d1 = np.load( 'debug/spatial_checkerboard_'+str( trial ) +\
                      '_full_'+str(1)+'.pdia.npy')
        pdfigs.append( (d0,d1,fig0, fig1) )
            

    return ( morans, pdfigs, patches )
    

if __name__ == "__main__":

    N = 5
    trials = 1
            
    binary = test_checkerboard_binary( N )

    # checker = test_checkerboard( N, 1 )
    # mats = checker[-1]
    # for i,m in enumerate( mats ):
    #     np.save( 'debug/spatial_checkerboard_patches_' +str(i), m )
