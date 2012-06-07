import inspect
import pylab
import numpy as np
import scipy

import csv

from sdt import SDT

def run_comparison(a, b, N,  pcolor):
    a_label = a.__doc__
    b_label = b.__doc__
##    N=10

    H,F,A,B = [],[],[],[]
    for f in xrange(N+1):
        H.append([])
        F.append([])
        A.append([])
        B.append([])

        for h in xrange(N+1):
            counts = {'HI':h, 'MI':N-h, 'FA':f, 'CR':N-f}
            
            H[-1].append(h)
            F[-1].append(f)
            A[-1].append(a(SDT(counts)))
            B[-1].append(b(SDT(counts)))
            
    H=np.array(H)
    F=np.array(F)
    A=np.array(A)
    B=np.array(B)

    A_vs_B  = np.array(B)
    A_vs_B -= np.min(A_vs_B)
    A_vs_B /= np.max(A_vs_B)
    A_vs_B -= A
            
    with open('%s-vs-%s,N=%i.csv'
              %(a_label,b_label,N), 'wb') as f:
        csv.writer(f).writerows(
            zip(*[H.flatten(),F.flatten(), A.flatten(),B.flatten()]))

    pylab.figure(figsize=(6,6))
    pylab.subplots_adjust(bottom=.15, top=.9, left=.15, right=.9)
    pylab.scatter(B,A, alpha=.5)
    pylab.ylabel(a_label)
    pylab.xlabel(b_label)
    pylab.savefig('%s-vs-%s,N=%i.png'%(a_label,b_label,N))
    pylab.close()

    if pcolor:
        pylab.figure(figsize=(4,8))
        pylab.subplots_adjust(bottom=.1, top=.95, hspace=.3)
        pylab.subplot(311, aspect='equal')
        pylab.title(a_label)
        pylab.pcolor(F,H,A)
        pylab.colorbar()

        pylab.subplot(312, aspect='equal')
        pylab.title(b_label)
        pylab.pcolor(F,H,B)
        pylab.colorbar()

        pylab.subplot(313, aspect='equal')
        pylab.title("norm(%s) - %s"%(b_label,a_label))
        pylab.pcolor(F,H,A_vs_B)
        pylab.colorbar()
        pylab.savefig('%s-vs-%s_pcolor,N=%i.png'
                      %(a_label,b_label,N),dpi=300)
        pylab.close()

if __name__ == '__main__':

    N=100
    for kwargs in [
        {'a':SDT.aprime,             'b':SDT.dprime,  'N':N, 'pcolor':True},
        {'a':SDT.amzs,               'b':SDT.dprime,  'N':N, 'pcolor':True},
        {'a':SDT.mutual_information, 'b':SDT.logbeta, 'N':N, 'pcolor':False},
        {'a':SDT.mutual_information, 'b':SDT.logbmz,  'N':N, 'pcolor':False},
        {'a':SDT.mutual_information, 'b':SDT.B,       'N':N, 'pcolor':False},
        {'a':SDT.mutual_information, 'b':SDT.c,       'N':N, 'pcolor':False},
        {'a':SDT.dprime,             'b':SDT.logbeta, 'N':N, 'pcolor':False},
        {'a':SDT.dprime,             'b':SDT.c,       'N':N, 'pcolor':False},
        {'a':SDT.amzs,               'b':SDT.logbmz,  'N':N, 'pcolor':False},
        {'a':SDT.aprime,             'b':SDT.bppd,    'N':N, 'pcolor':False},
        {'a':SDT.amzs,               'b':SDT.bppd,    'N':N, 'pcolor':False},
        {'a':SDT.bppd,               'b':SDT.c,       'N':N, 'pcolor':True},
        {'a':SDT.bppd,               'b':SDT.logbmz,  'N':N, 'pcolor':True},
        {'a':SDT.logbmz,             'b':SDT.logbeta, 'N':N, 'pcolor':True},
        {'a':SDT.logbmz,             'b':SDT.c,       'N':N, 'pcolor':True} ]:
        
        print kwargs
        run_comparison(**kwargs)
        
