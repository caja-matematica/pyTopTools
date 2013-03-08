from pyTopTools import chomp_betti as cb
import matplotlib.pyplot as plt
import pp


start = time.time()

stack_height = [10, 20]#

path = '/data/CT_Firn_Sample/output23-10-3/'
chomp_path = '/data/CT_Firn_Sample/chomp_files/'
prefix = 'K09b-23-10-'


chomp_path = '/sciclone/data10/jberwald/CT_Firn_Samples/chomp_files/'
path = '/sciclone/data10/jberwald/CT_Firn_Samples/output23-10-3/'
prefix = 'K09b-23-10-'

## Write Chomp-readable files for 3D blocks

#parallelize this stuff
ncpus = len( stack_height )
job_server = pp.Server( ncpus, ppservers=() )   
pool = []

bottom = 3200
top = 3220


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
        for base in [3310, 3320]:#range( 3200,  3400, height ):
            frames = []
            # list of frames to stack
            for x in range( height ):
                num = base + x
                frames.append( path + prefix + str( num ) + '.bmp' )

            print "    Stacking from base:", base
            stack = cb.stack_images( frames, height )
            cubfile = chomp_path + prefix[:-1] + \
                '_b' + str( base ) + \
                '_h' + str( height ) 

            # Convert bmp files to array, stack them, write them to
            # chomp-readable format.
            cb.array2chomp( stack, cubfile + '.cub' )

            # Now compute homology for each block
            cb.run_chomp( cubfile + '.cub', cubfile + '.hom'  )
            cb.extract_betti( cubfile + '.hom' )
        print ""

    print "Time:", time.time() - start

if 1:

    #for height in stack_height:
    height = 10
    dim = 1
    bettis = []
    for base in range( 3310,  3330, height ):
        betti_file = chomp_path + prefix[:-1] + \
            '_b' + str( base ) + \
            '_h' + str( height ) + \
            '.betti'
        print betti_file
        bnums = np.loadtxt( betti_file, dtype=np.uint8 )
        print bnums
        bettis.append( bnums[dim][1] )

    fig = plt.figure()
    ax = fig.gca()
    ax.plot( bettis, 'bo' )
    ax.set_xlabel( "Block number (height="+str(height)+")" )
    ax.set_ylabel( r"$\beta_{"+str( dim )+"}$" )
    ax.set_xlim( -1, len(bettis) )
    ax.set_ylim( min(bettis)-1, max(bettis)+1 )
    plt.show()
