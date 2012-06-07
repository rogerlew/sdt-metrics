import numpy as np
from __future__ import print_function

# Copyright (c) 2012, Roger Lew [see LICENSE.txt]

if __name__ == '__main__':
    from pyvttbl import DataFrame
    from sdt_metrics import SDT

    df = DataFrame()

    N = 10
    for HI in xrange(N+1):
        for FA in xrange(N+1):
            MI = N-HI
            CR = N-FA
            counts = {'HI':HI,'MI':MI,'CR':CR,'FA':FA}

            df.insert([('HI',HI),('FA',FA),
                       ('Accuracy',           SDT.accuracy(SDT(counts))),
                       ('PrecisionPPV',       SDT.precision(SDT(counts))),
                       ("Dprime",             SDT.dprime(SDT(counts))),
                       ("Aprime",             SDT.aprime(SDT(counts))),
                       ("Amzs",               SDT.amzs(SDT(counts))),
                       ('f1',                 SDT.f1(SDT(counts))),
                       ('MCC',                SDT.mcc(SDT(counts))),
                       ('Mutual_Information', SDT.mutual_info(SDT(counts))),
                       ('logbeta',            SDT.logbeta(SDT(counts))),
                       ('B',                  SDT.B(SDT(counts))),
                       ('C',                  SDT.c(SDT(counts))),
                       ("bpp",                SDT.bpp(SDT(counts))),
                       ("bppD",               SDT.bppd(SDT(counts))),
                       ("logbmz",             SDT.logbmz(SDT(counts)))])

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
    df.scatter_matrix(['HI','FA','Dprime','C','Amzs','B'],
                      diagonal='kde',trend='linear')
    
    
