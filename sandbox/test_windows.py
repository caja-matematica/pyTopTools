from pyTopTools import timeseries
import numpy as np


tmin = 4600
tmax = 5000
steps = 50
vrprefix = './data/genbif0_vr_window' + str(steps)
persname = './data/genbif0_pers_window' + str(steps)

a = np.arange( tmin, tmax, steps )
b = np.arange( tmin+steps, tmax+steps, steps )
endpts = zip( a,b )

windows = []
for s,t in endpts:
    windows.append( timeseries.Window( './genbif0.txt', s, t, 'vr' ) )

for window in windows:
    vrfile = vrprefix + '_' + str(window.tmin) +'.txt'
    window.convert2perseus( vrfile, stepsize=0.00005, nsteps=50 )
    window.compute_persistence( vrfile, persname )
    c = np.corrcoef( window.data[:-1], window.data[1:] )
    title = vrfile[:-4] + " -- corr( x(t), x(t+1) ) = " + str( c[0,1] )
    fig = window.draw_diagram( scale=True, title=title )
    fig.savefig( vrfile[:-4] + '.png' )
    
# Now, test near bifurcation
tmin = 8200
tmax = 8700
steps = 50
vrprefix = './data/genbif0_vr_window' + str(steps) +'.txt'
persname = './data/genbif0_pers_window' + str(steps)

# update the windows
a = np.arange( tmin, tmax, steps )
b = np.arange( tmin+steps, tmax+steps, steps )
endpts = zip( a,b )

windows = []
for s,t in endpts:
    windows.append( timeseries.Window( './genbif0.txt', s, t, 'vr' ) )

for window in windows:
    vrfile = vrprefix + '_' + str(window.tmin) +'.txt'
    window.convert2perseus( vrfile, stepsize=0.00005, nsteps=50  )
    window.compute_persistence( vrfile, persname )
    c = np.corrcoef( window.data[:-1], window.data[1:] )
    title = vrfile[:-4] + " -- corr(x(t), x(t+1)) = " + str( c[0,1] )
    fig = window.draw_diagram( scale=True, title=title )
    fig.savefig( vrfile[:-4] + '.png' )
    
