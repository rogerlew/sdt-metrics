from __future__ import print_function
from __future__ import division

# Copyright (c) 2012, Roger Lew [see LICENSE.txt]

import math

import pylab
import matplotlib
from matplotlib import pyplot as plt
import numpy as np
import scipy

import sdt_metrics

def metric_validation_plot(metric_name, levels=None, N=100, log=False):
    """
    produces a pcolor plot of the metric over ROC space

       arg:
          metric_name: string defining the metric to plot

       kwds:
          N: number of intervals for probability axes

          log: specifies whether log transform should be applied

    """
    func = getattr(sdt_metrics, metric_name)

    H,F,A,B = [],[],[],[]

    # pcolor thinks the H and F indices are bin edges so we need so
    # the shape of the final arrays need to be (N+2, N+2)
    for f in xrange(N+2):
        H.append([])
        F.append([])
        A.append([])
        for h in xrange(N+2):
            H[-1].append(h)
            F[-1].append(f)
            if h > N or f > N:
                # this would be on the top and right but get excluded
                # they are set to zero so it doesn't corrupt the colorbar
                A[-1].append(0)
            else:
                A[-1].append(func(h, N-h, N-f, f))

    F,H,A = np.array(F), np.array(H), np.array(A)

    ticks = np.linspace(0,N+1,5)
    matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
    f = plt.figure(figsize=(12,9))

    # sometimes matplotlib can be a pain! If you just use subplot the
    # contour plot is bigger than the pcolor plot because pcolor has a
    # colorbar and contour doesn't. To make it look right we have to
    # use gridspec.
    #
    # Nice technique to know in the long run.
    # http://matplotlib.sourceforge.net/users/gridspec.html
    gs = matplotlib.gridspec.GridSpec(1, 2, height_ratios=[1,1], width_ratios=[1,1])
    gs.update(wspace=0.01)
    ax1 = plt.subplot(gs[0], aspect='equal')
    if levels==None: # this makes it easier to figure out what the levels should be
        cs = pylab.contour(F[1:-1,1:-1],H[1:-1,1:-1],A[1:-1,1:-1],
                           colors='k')
    else:
        cs = pylab.contour(F[1:-1,1:-1],H[1:-1,1:-1],A[1:-1,1:-1],
                           levels=levels, colors='k')
        
    matplotlib.pyplot.clabel(cs, fontsize=8, inline=1)
    pylab.title(metric_name)
    pylab.xlim([0,N+1])
    pylab.ylim([0,N+1])
    pylab.xticks(ticks, ['%.2f'%(t/(N+1)) for t in ticks], rotation=30)
    pylab.yticks(ticks, ['%.2f'%(t/(N+1)) for t in ticks])
    pylab.ylabel('p(HI)')
    pylab.xlabel('p(FA)')
    
    
    ax2 = plt.subplot(gs[1], aspect='equal')
    pylab.title((metric_name, 'log(%s)'%metric_name)[log])
    pylab.pcolor(F,H,(A, np.log(A))[log])
    pylab.xlim([0,N+1])
    pylab.ylim([0,N+1])
    pylab.xticks(ticks, ['%.2f'%(t/(N+1)) for t in ticks], rotation=30)
    pylab.yticks(ticks, ['' for t in ticks])
    pylab.xlabel('p(FA)')
    pylab.colorbar()
    
    pylab.savefig('%s__lores.png'%metric_name,bbox_inches='tight',dpi=100)
    pylab.savefig('%s.png'%metric_name,bbox_inches='tight',dpi=300)
    pylab.savefig('%s.pdf'%metric_name,bbox_inches='tight')
    pylab.close()


def sequence(start, stop, step):
    return np.arange(start, stop + .5*step, step)

if __name__ == '__main__':
    metrics = [
               ('aprime',           sequence(.05,.95,.10)),
               ('amzs',             sequence(.05,.95,.10)),
               ('bpp',              sequence(-.9,.9,.2)),
               ('bph',              sequence(-.875,.875,.25)),
               ('bppd',             sequence(-.9,.9,.2)),
               ('bmz',              sequence(.1, 2.9, .4)),
               ('b',                sequence(.1,.9,.1)),
               ('dprime',           sequence(-3.,3.,.5)),
               ('beta',             sequence(.2, 2.9, .3)),
               ('c',                sequence(-1.5,1.5,.25)),
               ('accuracy',         sequence(.1,.9,.1)),
               ('mcc',              sequence(-.9,.9,.2)),
               ('precision',        sequence(.1,.9,.1)),
               ('recall',           sequence(.1,.9,.1)),
               ('f1',               sequence(.1,.9,.1)),
               ('ppv',              sequence(.1,.9,.1)),
               ('npv',              sequence(.1,.9,.1)),
               ('fdr',              sequence(.1,.9,.1)),
               ('sensitivity',      sequence(.1,.9,.1)),
               ('specificity',      sequence(.1,.9,.1)),
               ('mutual_info',      sequence(.0,.9,.1)),
               ('loglinear_bppd',   sequence(-.9,.9,.2)),
               ('loglinear_dprime', sequence(-3.,3.,.5)),
               ('loglinear_beta',   sequence(.2, 2.9, .3)),
               ('loglinear_c',      sequence(-1.5,1.5,.25)),
              ]

    for (metric,levels) in metrics:
        print(metric)
        log = metric in ['beta','bmz']
        metric_validation_plot(metric, levels, log=log)
