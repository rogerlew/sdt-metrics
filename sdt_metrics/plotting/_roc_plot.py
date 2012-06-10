from __future__ import print_function

# Copyright (c) 2012, Roger Lew [see LICENSE.txt]

import pylab
import numpy as np
import scipy

from numpy import pi
from scipy.stats import norm

import sdt_metrics
from .._sdt_metrics import ltqnorm,HI,MI,CR,FA
from ._mult_roc_plot import mult_roc_plot

_normdist = lambda x : np.exp(-x**2/2.)/np.sqrt(2*pi)

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
        
