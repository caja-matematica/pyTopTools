import persistence_tools as P
import perseus_wrap as pw


pt = P.Data( 'genbif0.txt', 'timeseries' )

block = 200
pname = '/Users/jberwald/github/local/caja-matematica/pyTopTools/sandbox/genbif_pers'
for t in range( 8000, 9000, 100 ):
    pt.convert2perseus( pname + str(t)+'.txt', scale=100, block=(t,t+block) )
    pw.perseus( pname + str(t)+'.txt', pname + str(t), dtype='cubtop' )
    fig = pw.plot_diagram( pname + str(t) + '_0.txt' )
    ax = fig.gca()
    ax.set_title( 'genbif -- t0=' +str( t ) +'; block='+str(block) )
    pw.plt.draw()
    fig.savefig( pname + str(t) + '_0.png' )
