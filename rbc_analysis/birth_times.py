import numpy
import os, re
import matplotlib.pyplot as plt
import cPickle as pkl
import rbc_histogram as rh

from collections import defaultdict

# def track_gens( genarr, err=40 ):
#     """
#     Given the stability of persistence diagrams, we expect that
#     different generators are no indistinghuishable whenever |b-b'| <
#     \eps and |d-d'| < \eps.

#     genarr : list of n x 2 numpy arrays consisting of all generators
#     with lifespan > noise (usually 40)

#     Latest robust gens located in /sciclone/data10/jberwald/RBC/cells/persout.
#     """
#     bt = []
#     for frame in genarr:
        
    

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
