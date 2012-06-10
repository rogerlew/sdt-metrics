from __future__ import print_function

# Copyright (c) 2012, Roger Lew [see LICENSE.txt]

import numpy as np
import math

from pyvttbl import DataFrame
from sdt_metrics import SDT, HI,MI,CR,FA, bppd, c, loglinear_c, loglinear_bppd

# something to put data in
df = DataFrame()

# make up data
P = 10
N = 10
for hi in xrange(P+1):
    for fa in xrange(N+1):
        df.insert([(HI,hi),(MI,N-hi),(CR,N-fa),(FA,fa)])

# calculate some metrics
df['bppd']=bppd(df[HI],df[MI],df[CR],df[FA])
df['loglinear_c']=loglinear_c(df[HI],df[MI],df[CR],df[FA])
df['loglinear_bppd']=loglinear_bppd(df[HI],df[MI],df[CR],df[FA])
df['c']=c(df[HI],df[MI],df[CR],df[FA])
df['phi']= 1.*df[HI]/(df[HI]+df[MI])
df['pllhi']= (df[HI] + .5)/(df[HI] + df[MI] + 1.)
                           
##    df.scatter_matrix(['HI','Accuracy','MCC',"Aprime",'Amzs',
##                       'f1',"Dprime",'Mutual_Information'],
##                      diagonal='kde',trend='linear')
##
##    df.scatter_matrix(['logbeta','logbmz','B','C','bpp','bppD'],
##                      diagonal='kde',trend='linear')
##    
##    df.scatter_matrix(['HI','FA','Mutual_Information','Dprime','C','logbeta'],
##                      diagonal='kde',trend='linear')
##
##    df.scatter_matrix(['HI','FA','Dprime','C','Amzs','B'],
##                      diagonal='kde',trend='linear')

df.scatter_matrix(['phi','pllhi'],
                  diagonal='kde',trend='linear')
                           
df.scatter_matrix(['c','loglinear_c','loglinear_bppd','bppd'],
                  diagonal='kde',trend='linear')

df2 = df.where('bppd != 1 and bppd != -1')

df2.scatter_matrix(['c','loglinear_c','loglinear_bppd','bppd'],
                   diagonal='kde',trend='linear',fname='filetered.png')
