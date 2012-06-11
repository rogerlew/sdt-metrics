from __future__ import print_function

# Copyright (c) 2012, Roger Lew [see LICENSE.txt]

import numpy as np
import math

from pyvttbl import DataFrame

import sdt_metrics
from sdt_metrics import SDT, HI,MI,CR,FA

# something to put data in
df = DataFrame()

# make up data
P = 10
N = 10
for hi in xrange(P+1):
    for fa in xrange(N+1):
        df.insert([(HI,hi),(MI,N-hi),(CR,N-fa),(FA,fa)])


    metrics = ['aprime',
               'amzs',
               'bpp',
               'bph',
               'bppd',
               'bmz',
               'b',
               'dprime',
               'beta',
               'c',
               'accuracy',
               'mcc',
               'precision',
               'recall',
               'f1',
               'ppv',
               'npv',
               'fdr',
               'sensitivity',
               'specificity',
               'mutual_info',
               'loglinear_bppd',
               'loglinear_dprime',
               'loglinear_beta',
               'loglinear_c']
    
for metric in metrics:
    func = getattr(sdt_metrics, metric)
    df[metric] = func(df[HI],df[MI],df[CR],df[FA])

df['log(beta)']=np.log(df['beta'])
df['log(bmz)']=np.log(df['bmz'])
df['sgn(bph)log(abs(bph)+1)']= np.log(np.abs(df['bph'])+1)
df['sgn(bph)log(abs(bph)+1)'][np.where(df['bph']<0)]*=-1
                           
##df.scatter_matrix(['dprime','aprime','amzs','accuracy'],
##                  diagonal='kde',trend='linear',alternate_labels=False)
##
##df.scatter_matrix(['mcc','f1','mutual_info','accuracy'],
##                  diagonal='kde',trend='linear',alternate_labels=False)
##
##df.scatter_matrix(['log(beta)','c','bppd','loglinear_bppd','log(bmz)'],
##                  diagonal='kde',trend='linear',alternate_labels=False)

# make pdfs
df.scatter_matrix(['dprime','aprime','amzs','accuracy'],
                  diagonal='kde',trend='linear',alternate_labels=False,
                  fname='scatter_matrix(dprime_X_aprime_X_amzs_X_accu'
                        'racy,diagonal=kde,alternate_labels=False).pdf')

df.scatter_matrix(['mcc','f1','mutual_info','accuracy'],
                  diagonal='kde',trend='linear',alternate_labels=False,
                  fname='scatter_matrix(mcc_X_f1_X_mutual_info_X_accu'
                        'racy,diagonal=kde,alternate_labels=False).pdf')

df.scatter_matrix(['log(beta)','c','bppd','loglinear_bppd','log(bmz)'],
                  diagonal='kde',trend='linear',alternate_labels=False,
                  fname='scatter_matrix(log(beta)_X_c_X_bppd_X_loglinear_'
                  'bppd_X_log(bmz),diagonal=kde,alternate_labels=False).pdf')
