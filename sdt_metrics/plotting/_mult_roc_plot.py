from __future__ import print_function

# Copyright (c) 2012, Roger Lew [see LICENSE.txt]

import pylab
import numpy as np
import scipy

from numpy import pi
from scipy.stats import norm

import sdt_metrics
from .._sdt_metrics import ltqnorm,HI,MI,CR,FA

_normdist = lambda x : np.exp(-x**2/2.)/np.sqrt(2*pi)

def mult_roc_plot(*args, **kwds):
    """
    Multiple Receiver Operating Characteristic (ROC) curvesPlot

       args:
          each arg should contain pairs of data and labels
          The data can be specified in 3 ways:
             1 argument:
                sdt_metrics.SDT object
          
             2 arguments:
                pHI
                pFA
              
             4 arguments:
                hit count
                miss count
                correction rejection count
                false alarm count

          labels should be strings and could contain latex

       kwds:
          metric: dprime, aprime, amzs (default is dprime)
          
          fname: outputname
          
          dpi: resolution of plot
    """
    #
    # bookkeeping
    #
    fname = kwds.get('fname','roc_plot.png')
    dpi = kwds.get('dpi',150)
    metric = kwds.get('metric','dprime')
    metric_func = getattr(sdt_metrics, metric)

    #
    # initialize the figure
    #
    pylab.figure(figsize=(5.5,5.5))
    pylab.subplots_adjust(left=.1, bottom=.1, top=.95, right=.95)
    pylab.plot([0,1],[0,1],'k:') # dotted diagonal line
    
    # used to change line plotting styles
    # gives 84 unique combinations. That should be enough, having more than that
    # would be pretty difficult to comprehend
    colors = 'bgrcmyk'
    linestyles = ['-','--','-.',':']
    markerstyles = 'hvs'
    
    #
    # loop through arguments and plot curves
    #
    # args should be (data, label) pairs    
    for j,(arg,label) in enumerate(args):
        # assume arg is an sdt object
        if isinstance(arg, sdt_metrics.SDT):
            sdt_obj = arg
            hi,mi,cr,fa = sdt_obj[HI],sdt_obj[MI],sdt_obj[CR],sdt_obj[FA]
            pHI,pFA = sdt_obj.p('HI'),sdt_obj.p('FA')
            metric_val = metric_func(hi,mi,cr,fa)
            
        # assume args are hit and false alarm rates
        elif len(arg) == 2:
            pHI,pFA = arg
            metric_val = metric_func(*arg)

        # assume args hit, miss, cr, and fa counts
        elif len(arg) == 4:
            hi,mi,cr,fa = arg
            sdt_obj = sdt_metrics.SDT(HI=hi,MI=mi,CR=cr,FA=fa)
            pHI,pFA = sdt_obj.p('HI'),sdt_obj.p('FA')
            metric_val = metric_func(*arg)

        #
        # build curve
        #
        if metric == 'dprime':
            # dprime we can handle quickly
            Z = np.linspace(-10,10,512)
            X = norm.cdf(Z-metric_val)
            Y = norm.cdf(Z)

        else:
            # this is probably not the best way to do this performance-wise...
            X = np.linspace(.001,.999,64).tolist()
            Y = []
            for pfa in X:
                func = lambda phi : abs(metric_func(phi,pfa)-metric_val)
                out = scipy.optimize.fminbound(func, x1=0., x2=1., xtol=1e-3)
                Y.append(out)

        #
        # plot data
        #
        
        # the actual roc curve
        pylab.plot(X, Y, 
                   c=colors[j%len(colors)],
                   ls=linestyles[j%len(linestyles)],
                   alpha=.6)
                    
        # the marker on the curve
        pylab.scatter([pFA], [pHI], 
                      c=colors[j%len(colors)],
                      marker=markerstyles[j%len(markerstyles)],
                      edgecolors='none',
                      s=40.,
                      alpha=.6)

        
        # a line with the color, linestyle, and marker style for 
        # just for the legend
        # (don't need to do this when there is only one arg supplied because
        # the legend doesn't show        
        if len(args) > 1:
            pylab.plot([-1,-2],[-1,-2],
                        alpha = .6,
                        c=colors[j%len(colors)],
                        ls=linestyles[j%len(linestyles)],
                        marker=markerstyles[j%len(markerstyles)],
                        markeredgewidth=0.,
                        label=label)

    #
    # do some final formatting
    #
    pylab.xlim([0,1])
    pylab.ylim([0,1])

    if len(args) > 1:
        pylab.legend(loc='lower right')
    #
    # save and close
    #
    pylab.savefig(fname, dpi=dpi)
    pylab.close()

def roc_plot(*args, **kwds):
    """
    Receiver Operating Characteristic (ROC) Plot

       args:
          1 argument:
             sdt_metrics.SDT object
          
          2 arguments:
             pHI
             pFA
             
          4 arguments:
             hit count
             miss count
             correction rejection count
             false alarm count

       kwds:
          metric: dprime, aprime, amzs (default is dprime)
          fname: outputname
          dpi: resolution of plot
    """
    # wrap mult_roc_plot
    if len(args) == 1:
        mult_roc_plot([args[0],''], **kwds)    
    else:
        mult_roc_plot([args,''], **kwds)
