from __future__ import print_function
from __future__ import division

# Copyright (c) 2012, Roger Lew [see LICENSE.txt]

# Python 2 to 3 workarounds
import sys
if sys.version_info[0] == 2:
    _strobj = basestring
    _xrange = xrange
elif sys.version_info[0] == 3:
    _strobj = str
    _xrange = range

if __name__ == '__main__':
    # requires pyvttbl which isn't explicitely listed as a dependency
    # so everything is going here to reduce likelihood of importation
    # issues

    #
    # import necessary libraries
    #
    import math

    import pylab
    import numpy as np

    from pyvttbl import DataFrame
    from sdt_metrics import SDT, HI,MI,CR,FA, bppd, c, loglinear_bppd

    #
    # Make up data
    #
    df = DataFrame() # more or less a dictionary with numpy arrays
    P = 10 # number of positive events
    N = 10 # number of negative events
    for hi in _xrange(P+1):
        for fa in _xrange(N+1):
            df.insert([(HI,hi),(MI,P-hi),(CR,N-fa),(FA,fa)])

    # calculate some metrics
    df['bppd'] = bppd(df[HI],df[MI],df[CR],df[FA])
    df['c'] = c(df[HI],df[MI],df[CR],df[FA])
    df['loglinear_bddp'] = loglinear_bppd(df[HI],df[MI],df[CR],df[FA])
    df['residuals'] = df['loglinear_bddp']-df['bppd']

    # exclude boundary conditions from this second table
    df2 = df.where('HI != 0 and MI != 0 and CR != 0 and FA != 0')

    #
    # perform analyses
    #
    df.scatter_matrix(['bppd','loglinear_bddp','c'], diagonal='kde')
    df.scatter_plot('c', 'bppd', trend='linear')
    df.scatter_plot('bppd', 'loglinear_bddp', trend='linear')
    df2.scatter_plot('bppd', 'loglinear_bddp', trend='linear',
                     fname='fscatter(bppd_X_bppdp,trend=linear,NoBoundaryCases).png')

    # plot of residuals
    pylab.figure(figsize=(8,4.5))
    pylab.subplots_adjust(right=.98)
    pylab.axhline(0,c='k',ls=':')
    pylab.scatter(df2['bppd'], df2['residuals'], alpha=.7)
    pylab.xlim([-1,1])
    pylab.savefig('loglinear_bddp_residuals.png')
    pylab.close()

    # plot of maximum value by N
    max_N = 500
    df3 = DataFrame()
    for n in _xrange(2,max_N+1):
        adj_p = .5/(n+1) # min value given number of positive or negative events
                         # assuming 50% prevalance
        max_val = bppd(adj_p, adj_p)
        df3.insert([('n', n), ('ll_bddp_max', max_val)])

    
    pylab.figure(figsize=(8,4.5))
    pylab.subplots_adjust(left=.10, bottom=.13, top=.98, right=.98)
    pylab.axhline(1,c='k',ls=':')
    pylab.plot(df3['n'],df3['ll_bddp_max'],alpha=.4)

    for n,v in zip(df3['n'],df3['ll_bddp_max']):
        # when n = 2 -> s = 1
        # when n = 500 -> s = .1111
        s = 1./math.log(n,2)
        pylab.scatter([n],[v], edgecolors='b', marker='x', s=s*80, alpha=s)

    # annotate asymptote
    for i,n in enumerate([2,3,5,10,20,30,50,100,200,300,500]):
            pylab.text(200, .985-i*.0053,
                       '%i  %.5f'%(n,df3['ll_bddp_max'][n-2]),
                       horizontalalignment='right')
    
    # format x-axis
    pylab.gca().set_xscale('log',basex=10)
    pylab.xlim([1.5,max_N])
    xticks = [2,3,5,10,20,30,50,100,200,300,500]
    pylab.xticks(xticks, [str(t) for t in xticks], rotation=30)
    pylab.xlabel('n = number positive events = number of negative events')

    # format y-axis
    pylab.ylim([.92-.005, 1.005])
    yticks = np.linspace(.92,1.,9).tolist()
    pylab.yticks(yticks, ['%.2f'%t for t in yticks])
    pylab.ylabel('maximum value of loglinear_bddp')
    
    pylab.savefig('ll_bddp_max(n).png',dpi=200)
    pylab.close()            

