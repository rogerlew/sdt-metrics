from __future__ import print_function

# Copyright (c) 2012, Roger Lew [see LICENSE.txt]
import math

import pylab
import numpy as np
import scipy

import sdt_metrics

def metric_validation_plot(metric_name, N=100, log=False, bph=False):
    """
    produces a pcolor plot of the metric over ROC space

       arg:
          metric_name: string defining the metric to plot

       kwds:
          N: number of intervals for probability axes

          log: specifies whether log transform should be applied

          bph: specifies the sgn(x)*log(abs(x)+1) transform
               used for bph
    """
    func = getattr(sdt_metrics, metric_name)
    if log:
        metric_name = 'log(%s)'%metric_name

    if bph:
        metric_name = 'sgn(%s)log(abs(%s)+1)'%(metric_name, metric_name)
        
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

    ticks = np.linspace(0,N+1,5)
    pylab.figure(figsize=(5,4))
    pylab.subplots_adjust(left=.2, bottom=.18, top=.93, right=.9)
    pylab.subplot(111, aspect='equal')
    pylab.title(metric_name)
    if bph: 
        mask = np.array([[(-1,1)[v>0] for v in L] for L in A])
        A = np.array([[math.log(abs(v)+1) for v in L] for L in A])
        pylab.pcolor(np.array(F), np.array(H), A*mask)
    elif log:
        pylab.pcolor(np.array(F), np.array(H), np.log(np.array(A)))
    else:
        pylab.pcolor(np.array(F), np.array(H), np.array(A))
    pylab.xlim([0,N+1])
    pylab.ylim([0,N+1])
    pylab.xticks(ticks, ['%.2f'%(t/(N+1)) for t in ticks], rotation=30)
    pylab.yticks(ticks, ['%.2f'%(t/(N+1)) for t in ticks])
    pylab.ylabel('p(HI)')
    pylab.xlabel('p(FA)')
    pylab.colorbar()
    pylab.savefig('%s.png'%metric_name)
    pylab.close()

if __name__ == '__main__':
    metrics = ['aprime',
               'amzs',
               'bpp',
               'bph',
##               'bppd',
##               'bmz',
##               'b',
##               'dprime',
##               'beta',
##               'c',
##               'accuracy',
##               'mcc',
##               'precision',
##               'recall',
##               'f1',
##               'ppv',
##               'npv',
##               'fdr',
##               'sensitivity',
##               'specificity',
##               'mutual_info',
##               'loglinear_bppd',
##               'loglinear_dprime',
##               'loglinear_beta',
               'loglinear_c']

    for metric in metrics:
        print(metric)
        log = metric in ['beta','bmz']
        bph = metric is 'bph'
        metric_validation_plot(metric, log=log, bph=bph)
