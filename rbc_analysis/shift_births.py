import birth_times as BT
import cPickle as pkl

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

with open( prefix + 'new_avg.pkl' ) as fh:
    newavg = pkl.load( fh )

with open( prefix + 'old_avg.pkl' ) as fh:
    oldavg = pkl.load( fh )

all_keys = oldavg.keys() + newavg.keys()

for cell_type in [old_name, new_name]:
    for val in lb:

        with open( prefix + cell_type + str(val) + '.pkl' ) as fh:
            cells = pkl.load( fh )
        ck = cells.keys()
        cell_vals = cells.values()

        cnames = [ s.split('/')[-2] for s in ck ]
        

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
                



        bt = BT.characteristic_birth_time( prefix + cell_type + str(val) + '.pkl',
                                           val,
                                           gen_num=gen_num )
    if avg_name in 
    with open( prefix + cell_type + str(val) + '_gen'+str(gen_num)+'.pkl', 'w' ) as fh:
        pkl.dump( bt, fh )
              

