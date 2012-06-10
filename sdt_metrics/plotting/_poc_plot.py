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

def poc_plot(*args, **kwds):
    # based on Ignacio Serrano-Pedraza Excel spreadsheet
    # "sdt_serranopedraza (version 1).xls"
    """
    Probability of Occurence Curves (POC)

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
          xmax: sets maximum x-limit
          
          fname: outputname
          
          dpi: resolution of plot
    """
    # process keyword arguments
    fname = kwds.get('fname','poc_plot.png')
    dpi = kwds.get('dpi',150)
    # xmax handled later (need dprime)
    
    # assume arg is an sdt object
    if len(args) == 1:
        sdt_obj = args[0]
        hi,mi,cr,fa = sdt_obj['HI'],sdt_obj['MI'],sdt_obj['CR'],sdt_obj['FA']
        dprime = sdt_obj.dprime()
        criterion = ltqnorm(1.-sdt_obj.p('HI')) + dprime

    # assume args are hit and false alarm rates
    elif len(args) == 2:
        pHI,pFA = args
        dprime = sdt_metrics.dprime(*args)
        criterion = ltqnorm(1.-pHI) + dprime
        
    # assume args hit, miss, cr, and fa counts
    elif len(args) == 4:
        hi,mi,cr,fa = args
        sdt_obj = sdt_metrics.SDT(HI=hi,MI=mi,CR=cr,FA=fa)
        dprime = sdt_obj.dprime()
        criterion = ltqnorm(1.-sdt_obj.p('HI')) + dprime

    # This represents the x-axis for the normal curves
    # in practice it is much longer than we need it, but we
    # want to be sure the end of the tail doesn't show with
    # high dprime values
    Z = np.linspace(-10,10,512)

    # this is a normal distribution -10 < Z < 10
    fxn = _normdist(Z)

    # initialize the figure (16/9 aspect ratio)
    pylab.figure(figsize=(8,4.5))
##    pylab.figure(figsize=(7,3))
    pylab.subplots_adjust(left=.08,right=.98)
    
    # plot the noise distribution
    pylab.plot(Z, fxn, 'r', lw=2., label = r'$f(x|n)$')#,alpha=.4)
    pylab.fill(Z, fxn, 'r', lw=0., alpha=.25)

    # plot the signal distribution    
    pylab.plot(Z + dprime, fxn, 'b--', lw=2., label = r'$f(x|s)$')#,alpha=.4)
    pylab.fill(Z + dprime, fxn, 'b', lw=0., alpha=.25)

    # annotate the peak of the signal distribution
    pylab.text(dprime, .41, r'$%.3f$'%dprime, horizontalalignment='center')

    # plot the criterion    
    pylab.axvline(criterion, color='g', ls=':', lw=2., label=r'$criterion$')

    # annotate the criterion 
    arrow_props = dict(facecolor='k', shrink=0.05, width=.25, headwidth=4., frac=.2)
    pylab.annotate(r'$%.3f$'%criterion, xy=(criterion, .44), xytext=(criterion-1.1, .46),
                   arrowprops=arrow_props, verticalalignment='bottom')

##    pylab.text(-2.5,.1, r'$sdt\_metrics$', fontsize=60) # used for Sphinx doc logo

    # format the y-axis
    pylab.ylim([0.0, 0.5])
    pylab.yticks([.0,.1,.2,.3,.4,.5])

    # add counts if available
    if len(args) != 2:
        pylab.plot([-3.6,-2.0],[.4,.4],'k')
        pylab.plot([-2.8,-2.8],[.35,.45],'k')
        pylab.text(-2.9, .41, '$%i$'%hi, horizontalalignment='right', verticalalignment='bottom')
        pylab.text(-2.9, .39, '$%i$'%fa, horizontalalignment='right', verticalalignment='top')
        pylab.text(-2.7, .41, '$%i$'%mi, horizontalalignment='left',  verticalalignment='bottom')
        pylab.text(-2.7, .39, '$%i$'%cr, horizontalalignment='left',  verticalalignment='top')

    # format the x-axis
    if 'xmax' in kwds:
        pylab.xlim([-4., kwds['xmax']])
    else:
        pylab.xlim([-4., 4.6+dprime])
        
    # show the legend
    # by default it is located in the upper right corner
    pylab.legend()

    # save and close figure
    pylab.savefig(fname,dpi=dpi)
    pylab.close()

