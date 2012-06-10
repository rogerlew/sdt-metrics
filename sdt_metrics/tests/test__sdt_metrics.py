# Copyright (c) 2011, Roger Lew [see LICENSE.txt]

"""
This unittest tests the sdt_metrics module.
"""

import sys
import unittest
import doctest
import random

from random import shuffle
from string import digits,ascii_lowercase

import sdt_metrics
from sdt_metrics import _S,SDT, HI,MI,FA,CR, mutual_info, aprime

class TestSDT__init__(unittest.TestCase):
    # Init test failure assertions
    def test0(self):
        with self.assertRaises(TypeError) as cm:
            SDT(42)

        self.assertEqual(str(cm.exception),
                 "'int' object is not iterable")
        
    def test1(self):
        with self.assertRaises(KeyError) as cm:
            SDT(one=1, two=2)

        self.assertEqual(str(cm.exception),"'two'")
        
    def test2(self):
        with self.assertRaises(KeyError) as cm:
            SDT([('one',1),('two',2)])

        self.assertEqual(str(cm.exception),"'one'")

    # test initialization signatures
    def test20(self):
        """SDT()"""
        self.assertEqual(repr(SDT()),"SDT()")

    def test21(self):
        """SDT(mapping)"""
        D = SDT(dict([(HI,10),(CR,9),(MI,1)]))
        self.assertEqual(repr(D), 'SDT(HI=10, MI=1, CR=9, FA=0)')

    def test22(self):
        """SDT(iterable)"""
        D = SDT([(HI,10),(CR,9),(MI,1)])
        self.assertEqual(repr(D), 'SDT(HI=10, MI=1, CR=9, FA=0)')

    def test23(self):
        """SDT(**kwargs)"""
        D = SDT(HI=10, MI=1, CR=9)
        self.assertEqual(repr(D), 'SDT(HI=10, MI=1, CR=9, FA=0)')

    def test24(self):
        """SDT(iterable)"""
        D = SDT([HI,HI,HI,FA,FA])
        self.assertEqual(repr(D), 'SDT(HI=3, MI=0, CR=0, FA=2)')
        
    def test25(self):
        """SDT(iterable, **kwargs), with overlapping key/values"""
        D = SDT(dict([(HI,10),(CR,9)]),MI=1)
        self.assertEqual(repr(D), 'SDT(HI=10, MI=1, CR=9, FA=0)')
        
    def test99(self):
        """Make sure that direct calls to update
           do not clear previous contents"""
        
        D = SDT(dict([(HI,10),(CR,9)]))
        D.__init__(MI=1)
        self.assertEqual(repr(D),'SDT(HI=10, MI=1, CR=9, FA=0)')

class TestSDT_clear(unittest.TestCase):
    def test0(self):        
        D = SDT(dict([(HI,10),(CR,9)]))
        D.clear()
        self.assertEqual(repr(D),'SDT()')
        
class TestSDT_delitem(unittest.TestCase):
    def test0(self):
        D = SDT(dict([(HI,10),(MI,1),(CR,9)]))
        del D[HI]
        self.assertEqual(repr(D),'SDT(HI=0, MI=1, CR=9, FA=0)')

class TestSDT_copy(unittest.TestCase):
    def test0(self):
        D = SDT(dict([(HI,10),(MI,1),(CR,9)]))
        D2 = D.copy()
        D2(HI)
        self.assertEqual(repr(D),'SDT(HI=10, MI=1, CR=9, FA=0)')

class TestSDT_fromkeys(unittest.TestCase):
    def test0(self):
        D = SDT(dict([(HI,10),(MI,1),(CR,9)]))
        
        with self.assertRaises(NotImplementedError) as cm:
            M=D.fromkeys([HI,HI])

        self.assertEqual(str(cm.exception),
             'SDT.fromkeys() is undefined.  Use SDT(iterable) instead.')
        
class TestSDT__setitem__(unittest.TestCase):
    def test0(self):
        D = SDT(HI=10,MI=1,CR=9)
        D[HI]=11
        self.assertEqual(repr(D),'SDT(HI=11, MI=1, CR=9, FA=0)')

    def test1(self):
        D = SDT(HI=10,MI=1,CR=9)
        with self.assertRaises(KeyError) as cm:
            D['AB']=11
        
        self.assertEqual(str(cm.exception),"'AB'")
        
class TestSDT_get(unittest.TestCase):
    def test0(self):
        self.assertEqual(SDT(HI=10,MI=1,CR=9).get(FA), 0) 
        
class TestSDT_setdefault(unittest.TestCase):
    def test0(self):
        D = SDT(dict([(HI,10),(MI,1),(CR,9)]))
        
        with self.assertRaises(NotImplementedError) as cm:
            M=D.setdefault('AB',3)

        self.assertEqual(str(cm.exception),
             'SDT.setdefault() is undefined.')
        
## update functions
class TestSDT_update(unittest.TestCase):
    
    def test0(self):
        L = SDT(HI=10,MI=1,CR=9)
        L.update(HI=11,FA=5)
        
        self.assertTrue(isinstance(L,SDT))        
        self.assertEqual(L,SDT(HI=21,MI=1,CR=9,FA=5)) # L is updated

class TestSDT_subtract(unittest.TestCase):
    # Init test failure assertions
    def test0(self):
        with self.assertRaises(TypeError) as cm:
            SDT().subtract(42)

        self.assertEqual(str(cm.exception),
                 "'int' object is not iterable")
        
    def test1(self):
        with self.assertRaises(KeyError) as cm:
            SDT().subtract(one=1, two=2)

        self.assertEqual(str(cm.exception),"'two'")
        
    def test2(self):
        with self.assertRaises(KeyError) as cm:
            SDT().subtract([('one',1),('two',2)])

        self.assertEqual(str(cm.exception),"'one'")

    # test initialization signatures
    def test20(self):
        """SDT()"""
        D = SDT()
        D.subtract()
        self.assertEqual(repr(D),"SDT()")

    def test21(self):
        """SDT(mapping)"""
        D = SDT()
        D.subtract(dict([(HI,10),(CR,9),(MI,1)]))
        self.assertEqual(repr(D), 'SDT(HI=-10, MI=-1, CR=-9, FA=0)')

    def test22(self):
        """SDT(iterable)"""
        D = SDT()
        D.subtract([(HI,10),(CR,9),(MI,1)])
        self.assertEqual(repr(D), 'SDT(HI=-10, MI=-1, CR=-9, FA=0)')

    def test23(self):
        """SDT(**kwargs)"""
        D = SDT()
        D.subtract(HI=10, MI=1, CR=9)
        self.assertEqual(repr(D), 'SDT(HI=-10, MI=-1, CR=-9, FA=0)')

    def test25(self):
        """SDT(iterable, **kwargs), with overlapping key/values"""
        D = SDT()
        D.subtract(dict([(HI,10),(CR,9)]),MI=1)
        self.assertEqual(repr(D), 'SDT(HI=-10, MI=-1, CR=-9, FA=0)')
                
    def test26(self):
        """SDT(iterable, **kwargs), with overlapping key/values"""
        D = SDT(dict([(HI,10),(CR,9)]),MI=1)
        D.update([HI,HI,HI,FA,FA])
        self.assertEqual(repr(D), 'SDT(HI=13, MI=1, CR=9, FA=2)')
        
class TestSDT__sub__(unittest.TestCase):
    def test0(self):
        L = SDT(HI=10,MI=1,CR=9)
        M = SDT(HI=11,FA=5)
        
        self.assertTrue(isinstance(L,SDT))        
        self.assertEqual(L - M, SDT(HI=0,MI=1,CR=9,FA=0))
        self.assertEqual(L, SDT(HI=10,MI=1,CR=9))
        self.assertEqual(M, SDT(HI=11,FA=5))
        
class TestSDT__or__(unittest.TestCase):
    def test0(self):
        L = SDT(HI=10,MI=1,CR=9)
        M = SDT(HI=11,FA=5)
        
        self.assertTrue(isinstance(L,SDT))        
        self.assertEqual(L | M, SDT(HI=11,MI=1,CR=9,FA=5))
        self.assertEqual(L, SDT(HI=10,MI=1,CR=9))
        self.assertEqual(M, SDT(HI=11,FA=5))
        
class TestSDT__and__(unittest.TestCase):
    def test0(self):
        L = SDT(HI=10,MI=1,CR=9)
        M = SDT(HI=11,FA=5)
        
        self.assertTrue(isinstance(L,SDT))        
        self.assertEqual(L & M, SDT(HI=10,MI=0,CR=0,FA=0))
        self.assertEqual(L, SDT(HI=10,MI=1,CR=9))
        self.assertEqual(M, SDT(HI=11,FA=5))

class TestSDT_items(unittest.TestCase):
    def test0(self):
        self.assertEqual(SDT().items(), zip([HI,MI,CR,FA],[0,0,0,0]))
        
class TestSDT_keys(unittest.TestCase):
    def test0(self):
        self.assertEqual(SDT().keys(), [HI,MI,CR,FA])
        
class TestSDT__iter__(unittest.TestCase):
    def test0(self):
        self.assertEqual(list(iter(SDT())), [HI,MI,CR,FA])
        
class TestSDT__call__(unittest.TestCase):
    def test0(self):
        D = SDT()
        D(HI)
        D(HI)
        D(HI)
        D(FA)
        D(FA)
        self.assertEqual(D, SDT(HI=3,FA=2))

    def test1(self):
        D = SDT()

        with self.assertRaises(KeyError) as cm:
            D('AB')

        self.assertEqual(str(cm.exception),"'AB'")

class TestSDT_count(unittest.TestCase):
    def test0(self):
        self.assertEqual(SDT().count(),0)

    def test1(self):
        D = SDT([HI,HI,HI])

        self.assertEqual(D.count(),3)

class TestSDT_p(unittest.TestCase):
    # http://www.linguistics.ucla.edu/faciliti/facilities/statistics/dprime.htm
    def test0(self):
        self.assertEqual(SDT(HI=20,MI=5,FA=10,CR=15).p(HI),0.8)

    def test1(self):
        self.assertEqual(SDT(HI=20,MI=5,FA=10,CR=15).p(MI),0.2)

    def test2(self):
        self.assertEqual(SDT(HI=20,MI=5,FA=10,CR=15).p(CR),0.6)

    def test3(self):
        self.assertEqual(SDT(HI=20,MI=5,FA=10,CR=15).p(FA),0.4)
        
class TestSDT_dprime(unittest.TestCase):
    # http://www.linguistics.ucla.edu/faciliti/facilities/statistics/dprime.htm
    def test0(self):
        d = SDT(HI=20,MI=5,FA=10,CR=15).dprime()
        self.assertAlmostEqual(d, 1.0949683355866173, 7)
        
class TestSDT_loglinear_dprime(unittest.TestCase):
    def test0(self):
        d = SDT(HI=20,MI=5,FA=10,CR=15).loglinear_dprime()
        self.assertAlmostEqual(d, 1.044498705934068, 7)

class TestSDT_c(unittest.TestCase):
    def test0(self):
        d = SDT(HI=20,MI=5,FA=10,CR=15).c()
        self.assertAlmostEqual(d, -0.29413706493331, 7)
        
class TestSDT_loglinear_c(unittest.TestCase):
    def test0(self):
        d = SDT(HI=20,MI=5,FA=10,CR=15).loglinear_c()
        self.assertAlmostEqual(d, -0.2788451754114444, 7)
        
class TestSDT_mutual_information(unittest.TestCase):
    # http://www.cs.ubc.ca/~murphyk/Teaching/CS340-Fall07/reading/rocHandout.pdf
    def test0(self):
        D = SDT(HI=80,MI=00,FA=10,CR=00)
        self.assertEqual(D.mutual_info(), 0)

    def test1(self):
        D = SDT(HI=80,MI=10,FA=00,CR=10)
        self.assertEqual(D.mutual_info(), 0.18645353727945904)

    def test2(self):
        D = SDT(HI=78,MI=12,FA=00,CR=10)
        self.assertEqual(D.mutual_info(), 0.17350094092658325)

class TestSDT_PPV(unittest.TestCase):
    # http://en.wikipedia.org/wiki/Positive_predictive_value
    def test0(self):
        D = SDT(HI=20,MI=180,FA=10,CR=1820)
        self.assertEqual(D.ppv(), 0.1)

class TestSDT_NPV(unittest.TestCase):
    # http://en.wikipedia.org/wiki/Positive_predictive_value
    def test0(self):
        D = SDT(HI=20,MI=180,FA=10,CR=1820)
        self.assertEqual(round(D.npv(),3), 0.995)

class TestSDT_specificity(unittest.TestCase):
    # http://en.wikipedia.org/wiki/Positive_predictive_value
    def test0(self):
        D = SDT(HI=20,MI=180,FA=10,CR=1820)
        self.assertEqual(round(D.specificity(),2), 0.91)

class TestSDT_sensitivity(unittest.TestCase):
    # http://en.wikipedia.org/wiki/Positive_predictive_value
    def test0(self):
        D = SDT(HI=20,MI=180,FA=10,CR=1820)
        self.assertEqual(round(D.sensitivity(),2), 0.67)

class Test_Singleton(unittest.TestCase):
    # test code to make sure sdt is really a singleton
    def test0(self):
        id1 = id(_S.getInstance().sdt)
        mutual_info.direct(43,54,65,34)
        id2 = id(_S.getInstance().sdt)
        mutual_info.direct(50,0,0,50)
        id3 = id(_S.getInstance().sdt)

        self.assertEqual(id1,id2)
        self.assertEqual(id2,id3)

class Test__vmethod_direct(unittest.TestCase):
    def test0(self):
        """float args"""
        self.assertEqual(aprime.direct(12,3,4,34),
                         SDT(HI=12,MI=3,CR=4,FA=34).aprime())
    def test1(self):
        """list args"""
        R = [SDT(HI=12,MI=3,CR=4,FA=34).aprime(),
             SDT(HI=12,MI=3,CR=4,FA=4).aprime()]

        D = aprime.direct([12,12],[3,3],[4,4],[34,4])

        for r,d in zip(R,D):
            self.assertAlmostEqual(r,d,7)

class Test__vmethod_prob(unittest.TestCase):
    def test0(self):
        """float args"""

        self.assertEqual(aprime.prob(12/15., 34/38.),
                         SDT(HI=12,MI=3,CR=4,FA=34).aprime())
    def test1(self):
        """list args"""
        R = [SDT(HI=12,MI=3,CR=4,FA=34).aprime(),
             SDT(HI=12,MI=3,CR=4,FA=4).aprime()]

        D = aprime.prob([12/15., 12/15.],
                        [34/38., 4/8.])

        for r,d in zip(R,D):
            self.assertAlmostEqual(r,d,7)

    def test2(self):
        """test _prob binding"""
        self.assertEqual(hasattr(mutual_info,'prob'), False)

class TestSDT__vmethod__call__(unittest.TestCase):
    def test1(self):
        self.assertEqual(str(aprime.prob([12/15., 12/15.], [34/38., 4/8.])),
                         str(aprime([12/15., 12/15.], [34/38., 4/8.])))

    def test2(self):
        self.assertEqual(str(aprime.prob(12/15., 34/38.)),
                         str(aprime(12/15., 34/38.)))            

    def test3(self):
        self.assertEqual(str(aprime.direct([12,12],[3,3],[4,4],[34,4])),
                         str(aprime([12,12],[3,3],[4,4],[34,4])))

    def test4(self):
        self.assertEqual(str(aprime.direct(12,3,4,34)),
                         str(aprime(12,3,4,34)))

class Test_plotting_poc_curve(unittest.TestCase):
    def test1(self):
        """given an SDT object"""
        sdt = SDT(HI=116, MI=30, CR=323, FA=80)
        sdt_metrics.plotting.poc_plot(sdt)
        
    def test2(self):
        """given probabilities"""
        sdt_metrics.plotting.poc_plot(.67, .43)
        
    def test3(self):
        """given an counts"""
        sdt_metrics.plotting.poc_plot(116, 30, 50, 50)

class Test_plotting_roc_curve(unittest.TestCase):
    def test1(self):
        """given an SDT object"""
        sdt = SDT(HI=116, MI=30, CR=323, FA=80)
        sdt_metrics.plotting.roc_plot(sdt)

    def test2(self):
        """given probabilities"""
        sdt_metrics.plotting.roc_plot(.67, .43)
        
    def test3(self):
        """given an counts"""
        sdt_metrics.plotting.roc_plot(116, 30, 50, 50)

class Test_plotting_mult_roc_curve(unittest.TestCase):
    def test1(self):
        """given an SDT object"""
        sdt_obj = SDT(HI=116, MI=30, CR=323, FA=80)
        sdt_probs = (.97,.22)
        sdt_counts = (76,67,80,65)
        sdt_metrics.plotting.mult_roc_plot((sdt_obj,  'from SDT object'),
                                           (sdt_probs, 'from probs'),
                                           (sdt_counts, 'from counts'),
                                           fname = 'mult_roc_example.png')
        

     
def suite():
    return unittest.TestSuite((
            unittest.makeSuite(TestSDT__init__),
            unittest.makeSuite(TestSDT_clear),
            unittest.makeSuite(TestSDT_delitem),
            unittest.makeSuite(TestSDT_get),
            unittest.makeSuite(TestSDT_setdefault),
            unittest.makeSuite(TestSDT_copy),
            unittest.makeSuite(TestSDT__setitem__),
            unittest.makeSuite(TestSDT_update),
            unittest.makeSuite(TestSDT_subtract),
            unittest.makeSuite(TestSDT__sub__),
            unittest.makeSuite(TestSDT__or__),
            unittest.makeSuite(TestSDT__and__),
            unittest.makeSuite(TestSDT_items),
            unittest.makeSuite(TestSDT_keys),
            unittest.makeSuite(TestSDT__iter__),
            unittest.makeSuite(TestSDT__call__),
            unittest.makeSuite(TestSDT_count),
            unittest.makeSuite(TestSDT_p),
            unittest.makeSuite(TestSDT_dprime),
            unittest.makeSuite(TestSDT_loglinear_dprime),
            unittest.makeSuite(TestSDT_c),
            unittest.makeSuite(TestSDT_loglinear_c),
            unittest.makeSuite(TestSDT_PPV),
            unittest.makeSuite(TestSDT_NPV),
            unittest.makeSuite(TestSDT_specificity),
            unittest.makeSuite(TestSDT_sensitivity),
            unittest.makeSuite(TestSDT_mutual_information),
            unittest.makeSuite(Test_Singleton),
            unittest.makeSuite(Test__vmethod_direct),
            unittest.makeSuite(Test__vmethod_prob),
            unittest.makeSuite(Test__vmethod_prob),
            unittest.makeSuite(Test_plotting_poc_curve),
            unittest.makeSuite(Test_plotting_roc_curve),
            unittest.makeSuite(Test_plotting_mult_roc_curve)
                              ))

if __name__ == "__main__":

    # run tests
    runner = unittest.TextTestRunner()
    runner.run(suite())
    
