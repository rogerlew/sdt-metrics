from __future__ import print_function
from __future__ import division

# Copyright (c) 2012, Roger Lew [see LICENSE.txt]

# Python 2 to 3 workarounds
import sys
if sys.version_info[0] == 2:
    _strobj = basestring
    _xrange = xrange
elif sys.version_info[0] == 3:
    _strobj = str
    _xrange = range
    
from _abcoll import Mapping
import math

from sdt_metrics.support.singletonmixin import Singleton

HI,MI,CR,FA = TP,TN,FN,FP = 'HI','MI','CR','FA'

"""
Acknowledgements
=================

  - SDT class scaffolded from Raymond Hettinger's Counter
    class in the collections module of the standard library

  - singletonmixin module is courtesy of Gary Robinson and is
    in public domain

  - ltqnorm was authored by Peter John Acklam

  - contributors to Python, Numpy, SciPy, and Matplotlib

Implementation Notes
====================

There are a couple of gotchas in the implementation that make this
code a little obfuscated. To help provide some insight I thought it
best to briefly discuss those gotchas here.

 Gotchas
 -------

 1. The inverse normal function used by dprime, beta, and c is
    undefined at 0 and 1 probabilities. The "standard correction" to
    set a probability of 0 to 1/(2N) and to set a probability of 1 to
    1-1/(2N). To make this as transparent to the user as possible I didn't
    want the function prototypes to be different for these metrics.
    This is part of the reason that the metrics are methods of the SDT
    class.

    When dprime.prob, beta.prob, and c.prob are passed 0 or 1
    probabilities an Exception is raised in _correction

    the loglinear_dprime, loglinear_beta, and loglinear_c methods
    use a different correction suggested by Hautus (1995).

 2. Some of the metrics (aprime, amzs, and bmz) have symmetry about the
    diagonal in ROC space where p(HI) == p(FA). Their implemenations are
    easier with recursion. These recurive metrics are implemented
    outside of the SDT class (_aprime, _amzs, _bmz, _bpp) and have
    wrappers inside the SDT class

 3. The SDT class was designed to be used as a score keeping
   datastructure. In some instances SDT metrics need to be calculated
   from frequency counts. In these instances having functions that
   take counts or probabilities makes more sense. To avoid implementing
   all the algorithms twice a singleton object (_S) with an SDT object
   (_S.sdt) is created when this class is imported. This prevents
   unnecessary instances of SDT being created.

   A factory class (_vmethod) let's us build vectorized versions of the
   metrics. To support taking probabilities directly the SDT class has
   a "toggle switch" (SDT._directmode) that changes how SDT.p() and
   SDT.count() behave.
   
      When self._directmode is True:
        p():     returns probabilities based on the frequency counts
        count(): returns the sum of the frequency counts

      When self._directmode is False
        p():     returns probabilities based on self.pHI and self.pFA.
        
                   -  These attributes are otherwise not used.
                   
                   -  They are set directly in _vmethod.prob
                 
        count(): returns None (to trigger exceptions for Gotcha 1)

  4. Some of the metrics need count data because the p(HI) and p(FA)
     do not convey the True prevalence rate. This is needed to calculate
     ppv, npv, f1, fdr, sensitivity, specificity, mutual_info, precision,
     and recall.

     To deal with this the _vmethod factory only adds the _prob function
     as a method of the metric if the metric supports it.

  5. Depending on the domain some users may want to use 'HI', 'FA',
     'CR', and 'FA' instead of 'TP', 'FN', and so on. This is why
     their is the global declaration:
         HI,MI,CR,FA = TP,FP,TN,FN = 'HI','MI','CR','FA'

         The idea is a user could do something like:
             from sdt import HI,MI,CR,FA

         or:
             from sdt import TP,FP,TN,FN

     and then use __call__ or __getitem__ functionality to update
     the counts with the terminology of their choice

     For Example:
         >>> from sdt import SDT, HI,MI,CR,FA
         >>> D = SDT()
         >>> D(HI)    # add a hit
         >>> D[HI]+=1 # adds another hit
"""

##
## SDT Support Functions
##

def _correction(v, N):
    """protects input to ltqnorm"""
    # used to protect input to ltqnorm
    # v is assumed to be a probability between 0 and 1
    if 0 < v < 1:
        return v
    elif N is None or v < 0 or v >1:
        raise ValueError('v should be >= 0 and <= 1')

    # at this point we know v must be 0 or 1 and N is not None
    return (1/(2*N), 1-1/(2*N))[int(v)]

def _isint(x):
    """returns True if x is an int, False otherwise"""
    try:
        int(x)
        return True
    except:
        return False

def ltqnorm( p ):
    # could be replaced with scipy.stats.norm.ppf,
    # but not including it makes it a pure python module
    #
    # as far as the licensing goes the algo. description page says:
    # "Use them as you wish, for whatever purpose."
    """
    Algorithm description at:
    http://home.online.no/~pjacklam/notes/invnorm/#The_distribution_function

    Function obtained from
    http://home.online.no/~pjacklam/notes/invnorm/impl/field/ltqnorm.txt
    
    Modified from the author's original perl code (original comments follow below)
    by dfield@yahoo-inc.com.  May 3, 2004.

    Lower tail quantile for standard normal distribution function.

    This function returns an approximation of the inverse cumulative
    standard normal distribution function.  I.e., given P, it returns
    an approximation to the X satisfying P = Pr{Z <= X} where Z is a
    random variable from the standard normal distribution.

    The algorithm uses a minimax approximation by rational functions
    and the result has a relative error whose absolute value is less
    than 1.15e-9.

    Author:      Peter John Acklam
    Time-stamp:  2000-07-19 18:26:14
    E-mail:      pjacklam@online.no
    WWW URL:     http://home.online.no/~pjacklam
    """

    if p <= 0 or p >= 1:
        # The original perl code exits here, we'll throw an exception instead
        raise ValueError( "Argument to ltqnorm %f must be in open interval (0,1)" % p )

    # Coefficients in rational approximations.
    a = (-3.969683028665376e+01,  2.209460984245205e+02, \
         -2.759285104469687e+02,  1.383577518672690e+02, \
         -3.066479806614716e+01,  2.506628277459239e+00)
    b = (-5.447609879822406e+01,  1.615858368580409e+02, \
         -1.556989798598866e+02,  6.680131188771972e+01, \
         -1.328068155288572e+01 )
    c = (-7.784894002430293e-03, -3.223964580411365e-01, \
         -2.400758277161838e+00, -2.549732539343734e+00, \
          4.374664141464968e+00,  2.938163982698783e+00)
    d = ( 7.784695709041462e-03,  3.224671290700398e-01, \
          2.445134137142996e+00,  3.754408661907416e+00)

    # Define break-points.
    plow  = 0.02425
    phigh = 1 - plow

    # Rational approximation for lower region:
    if p < plow:
       q  = math.sqrt(-2*math.log(p))
       return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / \
               ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)

    # Rational approximation for upper region:
    if phigh < p:
       q  = math.sqrt(-2*math.log(1-p))
       return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / \
                ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)

    # Rational approximation for central region:
    q = p - 0.5
    r = q*q
    return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q / \
           (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1)

def _aprime(pHI,pFA):
    """recursive private function for calculating A'"""
    pCR = 1 - pFA

    # use recursion to handle
    # cases below the diagonal defined by pHI == pFA
    if pFA > pHI:
        return 1 - _aprime(1-pHI ,1-pFA)

    # Pollack and Norman's (1964) A' measure
    # formula from Grier 1971
    if pHI == 0 or pFA == 1:
        # in both of these cases pHI == pFA
        return .5

    return .5 + (pHI - pFA)*(1 + pHI - pFA)/(4*pHI*(1 - pFA))

def _amzs(pHI,pFA):
    """recursive private function for calculating A_{MZS}"""

    # catch boundary cases
    if pHI == pFA == 0 or pHI == pFA == 1:
        return .5
    
    # use recursion to handle
    # cases below the diagonal defined by pHI == pFA
    if pFA > pHI:
        return 1 - _amzs(1-pHI, 1-pFA)
    
    # upper left quadrant
    # Mueller, Zhang (2006)
    if   pFA <= .5 <= pHI:
        return .75 + (pHI-pFA)/4 - pFA*(1-pHI)
    # cases above the diagonal defined by pHI == pFA
    # and not in the upper left quadrant
    elif pHI <= (1-pFA):
        if pHI == 0:
            return (3 + pHI - pFA)/4
        else:
            return (3 + pHI - pFA - pFA/pHI)/4
    else:
        if pFA == 1:            
            return (3 + pHI - pFA)/4
        else:
            return (3 + pHI - pFA - (1-pHI)/(1-pFA))/4

def _bmz(pHI,pFA):
    """recursive private function for calculating A_{MZS}"""
    
    # use recursion to handle
    # cases below the diagonal defined by pHI == pFA
    if pFA > pHI:
        return _bmz(1-pHI, 1-pFA)

    if   pFA <= .5 <= pHI:
        return (5-4*pHI)/(1+4*pFA)
    elif pFA < pHI < .5:
        return (pHI**2+pHI)/(pHI**2+pFA)
    elif .5 < pFA < pHI:
        return ((1-pFA)**2+(1-pHI))/((1-pFA)**2+(1-pFA))
    else: # pHI == pFA
        return 1.

def _bph(pHI,pFA):
    if pFA > pHI:
        return _bph(1-pHI, 1-pFA)
    
    if pHI <= 1 - pFA:
        num = pFA*(1-pFA)
        dem = pHI*(1-pHI)

        if dem == 0:
            return 1.
        else:
            return 1. - num/dem
    else:
        num = pHI*(1-pHI)
        dem = pFA*(1-pFA)

        if dem == 0:
            return -1.
        else:
            return num/dem - 1.
        
##
## SDT Class
##    
class SDT(dict):
    """
    A class for storing signal detection data

       Class is intended to be used a datastructure for keeping
       track of signal detection data. It also supports the
       the general functions of sdt_metrics.
    """
    # class is modelled from collections.Counter
    global HI,MI,CR,FA,TP,FP,TN,FN
    
    def __init__(self, iterable=None, **kwds):
        """Create a new, empty SDT object.  And if given, count elements
        from an input iterable.  Or, initialize the count from another mapping
        of elements to their counts.
        """
        super(SDT, self).__init__()
        self.update(iterable, **kwds)

        # used by the .prob methods
        self._directmode = True
        self.pHI = None
        self.pFA = None

    def keys(self):
        """returns list of event types"""
        return [HI,MI,CR,FA]

    def __iter__(self):
        """iterates over event types"""
        for k in [HI,MI,CR,FA]:
            yield k

    def items(self):
        """returns list of event type count pairs"""
        return [(k,self[k]) for k in self]
                
    def __missing__(self, key):
        """The count of elements not in the SDT is zero."""
        # Needed so that self[missing_item] does not raise KeyError
        return 0

    # Override dict methods where necessary
    @classmethod
    def fromkeys(cls, iterable, v=None):
        # There is no equivalent method for SDTs because setting v=1
        # means that no element can have a count greater than one.
        raise NotImplementedError(
            'SDT.fromkeys() is undefined.  Use SDT(iterable) instead.')

    def setdefault(self, key, value):
        raise NotImplementedError('SDT.setdefault() is undefined.')
    
    def __call__(self, event):
        """adds event to object

           D = SDT()
           ``D(HI) <==> D.update([HI]) <==> D[HI]+=1``
        """
        if event not in self.keys():
            raise KeyError(event)
        self[event] += 1
        
    def update(self, iterable=None, **kwds):
        """Like dict.update() but add counts instead of replacing them.

        Source can be an iterable, a dictionary, or another SDT instance.
        """
        # The regular dict.update() operation makes no sense here because the
        # replace behavior results in the some of original untouched counts
        # being mixed-in with all of the other counts for a mismash that
        # doesn't have a straight-forward interpretation in most counting
        # contexts.  Instead, we implement straight-addition.  Both the inputs
        # and outputs are allowed to contain zero and negative counts.

        if iterable is not None:
            if hasattr(iterable, '__getitem__'):
                if hasattr(iterable, 'keys'):
                    for elem, count in iterable.iteritems():
                        if elem not in self.keys():
                            raise KeyError(elem)
                        
                        self[elem] = self.get(elem, 0) + count
                else:
                    for val in iterable:
                        if isinstance(val, str):
                            if val not in self.keys():
                                raise KeyError(val)
                            self[val] = self.get(val, 0) + 1
                        else:
                            elem, count = val
                            if elem not in self.keys():
                                raise KeyError(elem)
                            self[elem] = self.get(elem, 0) + count
            else:
                raise TypeError("'%s' object is not iterable"
                                % type(iterable).__name__)
        if kwds:
            self.update(kwds)
            
        for elem in self:
            if not self.has_key(elem):
                self[elem] = 0

    def __setitem__(self, key, value):
        """``sdt.__setitem(key) <==> sdt[key]``"""
        if key not in self.keys():
            raise KeyError(key)
        else:
            super(SDT, self).__setitem__(key, value)

    def subtract(self, iterable=None, **kwds):
        """
        Like dict.update() but subtracts counts instead of replacing them.
        Counts can be reduced below zero.  Both the inputs and outputs are
        allowed to contain zero and negative counts.

        Source can be an iterable, a dictionary, or another SDT instance.
        """
        if iterable is not None:
            if hasattr(iterable, '__getitem__'):
                if hasattr(iterable, 'keys'):
                    for elem, count in iterable.iteritems():
                        if elem not in self.keys():
                            raise KeyError(elem)
                        
                        self[elem] = self.get(elem, 0) - count
                else:
                    for val in iterable:
                        if isinstance(val, str):
                            self[val] = self.get(val, 0) - 1
                        else:
                            elem, count = val
                            if elem not in self.keys():
                                raise KeyError(elem)
                            self[elem] = self.get(elem, 0) - count
            else:
                raise TypeError("'%s' object is not iterable"
                                % type(iterable).__name__)
        if kwds:
            self.subtract(kwds)

    def copy(self):
        """Return a shallow copy."""
        return self.__class__(self)

    def __reduce__(self):
        return self.__class__, (dict(self),)

    def __delitem__(self, elem):
        """
        Like dict.__delitem__() but does not raise KeyError for missing values.
        """
        if elem in self:
            super(SDT, self).__delitem__(elem)

    def __repr__(self):
        if self.count()==0:
            return '%s()' % self.__class__.__name__
        items = ', '.join(['%s=%i'%(k,self[k]) for k in self])
        return '%s(%s)' % (self.__class__.__name__, items)

    # Multiset-style mathematical operations discussed in:
    #       Knuth TAOCP Volume II section 4.6.3 exercise 19
    #       and at http://en.wikipedia.org/wiki/Multiset
    #
    # Outputs guaranteed to only include positive counts.
    #
    # To strip negative and zero counts, add-in an empty SDT:
    #       c += SDT()

    def __add__(self, other):
        """Add counts from two SDTs."""
        if not isinstance(other, SDT):
            return NotImplemented
        result = SDT()
        for elem, count in self.items():
            newcount = count + other[elem]
            if newcount > 0:
                result[elem] = newcount
        for elem, count in other.items():
            if elem not in self and count > 0:
                result[elem] = count
        return result

    def __sub__(self, other):
        """Subtract count, but keep only results with positive counts."""
        if not isinstance(other, SDT):
            return NotImplemented

        diffs = []
        for k in self:
            d = self[k]-other[k]
            if d < 0:
                diffs.append(0)
            else:
                diffs.append(d)
                    
        return SDT(zip(self.keys(), diffs))

    def __or__(self, other): # overloads |
        """Union is the maximum of value in either of the input SDTs."""
        if not isinstance(other, SDT):
            return NotImplemented
        
        return SDT([(k,max(self[k],other[k])) for k in self])

    def __and__(self, other): # overloads &
        """Intersection is the minimum of corresponding counts."""
        if not isinstance(other, SDT):
            return NotImplemented
        
        return SDT([(k,min(self[k],other[k])) for k in self])

    def count(self):
        """returns count of events"""
        if self._directmode:
            return sum(v for v in self.values())
        else:
            return None
    
    def p(self, elem):
        """returns probability of event type"""
        if elem not in self.keys():
            raise KeyError(elem)

        if self._directmode:
            if   elem == HI : return self[HI]/float(self[HI] + self[MI])
            elif elem == MI : return self[MI]/float(self[HI] + self[MI])
            elif elem == CR : return self[CR]/float(self[CR] + self[FA])
            else            : return self[FA]/float(self[CR] + self[FA])

        # this branch is only used by the _vmethod.probs
        else: 
            if   elem == HI : return self.pHI
            elif elem == MI : return 1-self.pHI
            elif elem == CR : return 1-self.pFA
            else            : return self.pFA
            
    def aprime(self):
        """
        A': Non-parametric measure of sensitivity

          Devised by Pollack and Norman [1]_ but was reformalated and
          popularized by Grier [2]_. Smith [3]_ provides a nice
          historical overview, literature review, and analysis of A'.

          .. seealso:: :py:meth:`amzs`

          .. [1] Pollack, I., Norman, D. A. (1964). A non-parametric analysis
                 of recognition experiemnts. Psychonomic Sicence 1, 125-126.

          .. [2] Grier, J. B. (1971). Nonparametric indexes for sensitivity and
                 bias: Computing formulas. Psychological Bulletin, 75, 424-429.

          .. [3] Smith, W. D. (1995). Clarification of Sensitivity Measure A'.
                 Journal of Mathematical Psychology 39, 82-89.        
        """
        return _aprime(self.p(HI), self.p(FA))
    
    def amzs(self):
        """
        A: Zhang and Mueller's ROC-Based Measure of Sensitivity

          Smith [1]_ remedied a common confusion with A' a suggested an
          improved measure A'' (not implemented). Zhang and Mueller [2]_ found
          that Smith had a mathematical error and properly formulated a new
          nonparametric measure of sensitivity. They called this measure A.
          Here it is called :py:meth:`amzs` in honor of Mueller, Zhang and Smith.        

          .. seealso:: :py:meth:`aprime`
          
          .. [1] Zhang, J., and Mueller, S. T. (2005). A note on roc analysis
                 and non-parametric estimate of sensitivity. Psychometrika, 70,
                 145-154.
          
          .. [2] Smith, W. D. (1995). Clarification of Sensitivity Measure A'.
                 Journal of Mathematical Psychology 39, 82-89.
        """
        return _amzs(self.p(HI),self.p(FA))

    def bpp(self):
        """
        b'': Grier [1]_ nonparametric measure of response bias

          .. seealso::
             :py:meth:`bppd`
             :py:meth:`loglinear_bppd`
          
          .. [1] Grier, J. B. (1971). Nonparametric indexes for sensitivity and
                 bias: Computing formulas. Psychological Bulletin, 75, 424-429.
        """
        pHI,pFA = self.p(HI),self.p(FA)
        if pHI >= pFA:
            num = pHI*(1-pHI) - pFA*(1-pFA)
            dem = pHI*(1-pHI) + pFA*(1-pFA)
        else:
            num = pFA*(1-pFA) - pHI*(1-pHI)
            dem = pFA*(1-pFA) + pHI*(1-pHI)

        if dem == 0:
            return 0.
        else:
            return num/dem
        
    def bph(self):
        """
        B'H: nonparametric measure of response bias

          Developed by Hodos [1]_. See, Warm, Dember, and Howe, [2]_ compared
          several nonparmetric measures of response bias and found bppd to
          be empirically better.

          .. seealso:: :py:meth:`bppd`
          
          .. [1] Hodos, W. (1970). A nonparametric index of response bias for
                 use in detection and recognition experiments. Psychological
                 Bulletin, 74, 351-354.

          .. [2] See, J. E., Warm, J. S., Dember, W. N., and Howe, S. R. (1997).
                 Vigilance and signal detection theory: An empirical evaluation
                 of five measures of response bias.
        """
        pHI,pFA = self.p(HI),self.p(FA)
        return _bph(pHI,pFA)

    def bppd(self):
        """
        B''d: nonparametric measure of response bias

          First developed by Donaldson [1]_. See, Warm, Dember, and
          Howe, [2]_ compared several nonparmetric measures of
          response bias and endorse this measure when the assumption of a
          gaussian noise distribution does not hold.

          .. seealso:: :py:meth:`loglinear_bppd`
          
          .. [1] Donaldson, W. (1992). Measuring recognition memory. Journal of
                 Experimental Psychology: General, 121, 275277.

          .. [2] See, J. E., Warm, J. S., Dember, W. N., and Howe, S. R. (1997).
                 Vigilance and signal detection theory: An empirical evaluation
                 of five measures of response bias.
        """
        pHI,pFA = self.p(HI),self.p(FA)
        num = ((1.-pHI)*(1.-pFA)-pHI*pFA)
        dem = ((1.-pHI)*(1.-pFA)+pHI*pFA)
        if dem == 0:
            return 0.
            
        return num / dem

    def loglinear_bppd(self):
        """
        loglinear B''d: nonparametric measure of response bias

          Corrects boundary problems with bppd. by applying Hautus's [1]_ loglinear
          transformation to counts before calculating bppd. For more information see
          the :doc:`loglinear_bppd_analysis`
          
          .. note:: This metric cannot be calculated from probabilities

          .. seealso:: :py:meth:`bppd`
              
          .. [1] Hautus, M. (1995). Corrections for extreme proportions and
                 their biasing effects on estimated values of d'. Behavior
                 Research Methods, Instruments, and Computers, 27, 46-51.

        """
        pHI = (self[HI] + 0.5)/(self[HI] + self[MI] + 1)
        pFA = (self[FA] + 0.5)/(self[CR] + self[FA] + 1)
        
        num = ((1.-pHI)*(1.-pFA)-pHI*pFA)
        dem = ((1.-pHI)*(1.-pFA)+pHI*pFA)

        return num / dem
        
    def bmz(self):
        """
        b: Zhang and Mueller's [1]_ measure of decision bias
        
          .. [1] Zhang, J., and Mueller, S. T. (2005). A note on roc analysis
                 and non-parametric estimate of sensitivity. Psychometrika, 70,
                 145-154.
        """
        return _bmz(self.p(HI),self.p(FA))

    def b(self):
        """
        B: 0.5*p(HI) + 0.5p(FA)

          Bar-Napkin measure of Response Bias.
          Implemented for comparitive purposes.
        """
        return 0.5*self.p(HI) + 0.5*self.p(FA)
    
    def dprime(self):
        """
        d': parametric measure of sensitivity

          Extremely popular measure adapted by from communication
          engineering by psychologists in the 1950s [1]_. Most notable text
          is Green and Swets  [2]_. Calculation uses the formula given
          by Macmillan [3]_. Extreme probabilities of 0 and 1 are
          treated using the correction suggested by Macmillan and Kaplan
          [4]_. Rates of 0 are replaced with 1/(2n) and rates of 1 are
          replaced wtih 1 - 1/(2n). This approach has been shown to be biased
          (Miller [5]_). Hautus's [6]_ loglinear approach may
          be preferred. For a more extensive overview on treating extreme values
          see Stanislaw and Todorov [7]_. 

          .. seealso:: :py:meth:`loglinear_dprime`
          
          .. [1] Szalma, J. L., and Hancock, P. A. Signal detection theory. Class
                 Lecture Notes. http://bit.ly/KIyKkt

          .. [2] Green, D. M., and Swets J. A. (1996/1988). Signal Detection
                 theory and psychophysics, reprint edition. Los Altos, CA:
                 Penisula Publishing.

          .. [3] Macmillan, N. A. (1993). Signal detection theory as data analysis
                 method and psychological decision model. In G. Keren & C. Lewis
                 (Eds.), A handbook for data analysis in the behavioral sciences:
                 Methodological issues (pp. 21-57). Hillsdale, NJ: Erlbaum.

          .. [4] Macmillan, N. A., and Kaplan, H. L. (1985). Detection theory
                 analysis of group data: Estimating sensitivity from average hit
                 and false-alarm rates. Psychological Bulletin, 98, 185-199.

          .. [5] Miller, J. (1996). The sampling distribution of d'. Perception
                 and Psychophysics, 58, 65-72.
              
          .. [6] Hautus, M. (1995). Corrections for extreme proportions and
                 their biasing effects on estimated values of d'. Behavior
                 Research Methods, Instruments, and Computers, 27, 46-51.

          .. [7] Stanislaw H. and Todorov N. (1999). Calculation of signal
                 detection theory measures. Behavorial Research Methods,
                 Instruments, and Computers, 31 (1), 137-149.
        """
        N = self.count()
        return ltqnorm(_correction(self.p(HI),N)) - \
               ltqnorm(_correction(self.p(FA),N))

    def loglinear_dprime(self):
        """
        loglinear d': parametric measure of sensitivity

          Applies Hautus's [1]_ loglinear transformation to counts
          before calculating dprime. (Transformation is applied to all
          values)

          .. note:: This metric cannot be calculated from probabilities

          .. seealso:: :py:meth:`dprime`
              
          .. [1] Hautus, M. (1995). Corrections for extreme proportions and
                 their biasing effects on estimated values of d'. Behavior
                 Research Methods, Instruments, and Computers, 27, 46-51.
        """
        pHI = (self[HI] + 0.5)/(self[HI] + self[MI] + 1)
        pFA = (self[FA] + 0.5)/(self[CR] + self[FA] + 1)
        
        return ltqnorm(pHI) - ltqnorm(pFA)
    
    
    def beta(self):
        """
        beta: classic parametric measure of response bias [1]_.

          Use of :py:meth:`c` (or :py:meth:`loglinear_c`) as a parametric
          measure of response bias may be preferred.

          This funciton applies
          Macmillan and Kaplan [2]_ correction to extreme values. 

          .. seealso::
             :py:meth:`c`
             :py:meth:`loglinear_c`
             :py:meth:`loglinear_beta`
        
          .. [1] Green, D. M., and Swets J. A. (1996/1988). Signal Detection
                 theory and psychophysics, reprint edition. Los Altos, CA:
                 Penisula Publishing.
              
          .. [2] Macmillan, N. A., and Kaplan, H. L. (1985). Detection theory
                 analysis of group data: Estimating sensitivity from average hit
                 and false-alarm rates. Psychological Bulletin, 98, 185-199.
        """
        N = self.count()
        zhr = ltqnorm(_correction(self.p(HI),N))
        zfar = ltqnorm(_correction(self.p(FA),N))
        return math.exp(-zhr*zhr/2 + zfar*zfar/2)


    def loglinear_beta(self):
        """
        beta: classic parametric measure of response bias [1]_.

          Use of :py:meth:`c` (or :py:meth:`loglinear_c`) as a parametric
          measure of response bias may be preferred.

          Applies Hautus's [2]_ loglinear transformation to counts before
          calculating beta. (Transformation is applied to all values)

          .. note:: This metric cannot be calculated from probabilities

          .. seealso::
             :py:meth:`c`
             :py:meth:`loglinear_c`
             :py:meth:`beta`


          .. [1] Green, D. M., and Swets J. A. (1996/1988). Signal Detection
                 theory and psychophysics, reprint edition. Los Altos, CA:
                 Penisula Publishing.
              
          .. [2] Hautus, M. (1995). Corrections for extreme proportions and
                 their biasing effects on estimated values of d'. Behavior
                 Research Methods, Instruments, and Computers, 27, 46-51.
        """
        pHI = (self[HI] + 0.5)/(self[HI] + self[MI] + 1)
        pFA = (self[FA] + 0.5)/(self[CR] + self[FA] + 1)
        zhr = ltqnorm(pHI)
        zfar = ltqnorm(pFA)
        return math.exp(-zhr*zhr/2 + zfar*zfar/2)
        
    def c(self):
        """
        c: parametric measure of response bias

          Generally recommended as a better measure than beta [1]_, [2]_,
          [3]_. First reason being that d' and c are independent [4]_.
          Forumula from Macmillan  [5]_.  Applies Macmillan and
          Kaplan correction to extreme values. [6]_.

          .. seealso::
             :py:meth:`dprime`
             :py:meth:`loglinear_c`

        
          .. [1] Banks W. P. (1970). Signal detection theory and human memory.
                 Psychological Bulletin, 74, 81-99.

          .. [2] Macmillan, N. A., and Creelman, C. D. (1990). Response bias:
                 Characteristics of detection theory, threshold theory, and
                 nonparametric indexes. Psychological Bulletin, 107, 401-413.

          .. [3] Snodgrass, J. G., and Corwin, J. (1988). Pragmatics of
                 measuring recognition memory: Applications to dementia and
                 amnesia. Journal of Experimental Psychology: General, 117, 34-50.

          .. [4] Ingham, J. G. (1970). Individual differences in signal detection.
                 Acta Psychologica, 34, 39-50.
          
          .. [5] Macmillan, N. A. (1993). Signal detection theory as data analysis
                 method and psychological decision model. In G. Keren and C. Lewis
                 (Eds.), A handbook for data analysis in the behavioral sciences:
                 Methodological issues (pp. 21-57). Hillsdale, NJ: Erlbaum.
              
          .. [6] Macmillan, N. A., and Kaplan, H. L. (1985). Detection theory
                 analysis of group data: Estimating sensitivity from average hit
                 and false-alarm rates. Psychological Bulletin, 98, 185-199.
        """
        N = self.count()
        return -1.*(.5*ltqnorm(_correction(self.p(HI),N)) + \
                    .5*ltqnorm(_correction(self.p(FA),N)))

    def loglinear_c(self):
        """
        c: parametric measure of response bias

          Applies Hautus's [1]_ loglinear transformation to counts
          before calculating beta. (Transformation is applied to all
          values)
          
          .. note:: This metric cannot be calculated from probabilities

          .. seealso::
             :py:meth:`c`
              
          .. [1] Hautus, M. (1995). Corrections for extreme proportions and
                 their biasing effects on estimated values of d'. Behavior
                 Research Methods, Instruments, and Computers, 27, 46-51.
        """
        pHI = (self[HI] + 0.5)/(self[HI] + self[MI] + 1)
        pFA = (self[FA] + 0.5)/(self[CR] + self[FA] + 1)
        
        return -1.*(.5*ltqnorm(pHI) + .5*ltqnorm(pFA))
    
    def accuracy(self):
        """
        accuracy: (1.+p(HI) - p(FA)) / 2.
        """
        return (1.+self.p(HI)-self.p(FA))/2.
    
    def mcc(self):
        """
        Matthews [1]_ correlation coefficient

          .. [1] Matthews, B.W. (1975). Comparison of the predicted and observed
                 secondary structure of T4 phage lysozyme. Biochim. Biophys. Acta,
                 405, 442-451.
        """
        pHI,pFA,pCR,pMI = self.p(HI),self.p(FA),self.p(CR),self.p(MI)
        num = pHI * pCR - pFA * pMI
        dem = math.sqrt((pHI + pFA)*( pHI + pMI )*( pCR + pFA )*( pCR + pMI ))
        
        if dem == 0:
            return 0
        return num / dem
    
    def precision(self):
        """
        precision
        
          .. note:: This metric cannot be calculated from probabilities

          ``sdt.precision() <==> sdt.ppv()``
        """
        return self.ppv()
    
    def recall(self):
        """
        recall
        
          .. note:: This metric cannot be calculated from probabilities

          ``sdt.recall() <==> sdt.sensitivity()``
        """
        return self.sensitivity()
    
    def f1(self):
        """
        f1: a.k.a. F1-Score, F-measure

          The harmonic mean of precision and recall [1]_
          
          .. note:: This metric cannot be calculated from probabilities

          .. [1] van Rijsbergen, C. J. (1979). Information Retrieval
                 (2nd ed.). Butterworth. http://bit.ly/8yAjR
        """
        precision,recall = self.precision(),self.recall()
        num = (2. * precision * recall)
        dem = (precision + recall)
        if dem == 0:
            return 0
        return num/dem
            
    def ppv(self):
        """
        positive predictive value: TP / (TP + FP)
        
          .. note:: This metric cannot be calculated from probabilities
        
          ``sdt.precision() <==> sdt.ppv()``
        """
        num = self[TP]
        dem = self[TP] + self[FP]
        if dem == 0:
            return 0.
        return num / float(dem)
    
    def npv(self):
        """
        negative predictive value: TN / (TN + FN)
        
          .. note:: This metric cannot be calculated from probabilities
        """
        num = self[TN]
        dem = self[TN] + self[FN]
        if dem == 0:
            return 0
        return num / float(dem)
    
    def fdr(self):
        """
        false discovery rate: FP / (TP + FP)
        
          .. note:: This metric cannot be calculated from probabilities
        """
        num = self[FP]
        dem = self[FP] + self[TP]
        if dem == 0:
            return 0
        return num / float(dem)

    def sensitivity(self):
        """
        sensitivity: TP / (TP + FN)
        
          .. note:: This metric cannot be calculated from probabilities
        
          ``sdt.recall() <==> sdt.sensitivity()``
        """
        num = self[TP]
        dem = self[TP] + self[FN]
        if dem == 0:
            return 0
        return num / float(dem)

    def specificity(self):
        """
        specificity: TN / (TN + FP)
        
          .. note:: This metric cannot be calculated from probabilities
        """
        num = self[TN]
        dem = self[TN] + self[FP]
        if dem == 0:
            return 0
        return num / float(dem)
        
    def mutual_info(self):
        # http://www.cs.ubc.ca/~murphyk/Teaching/CS340-Fall07/reading/rocHandout.pdf
        """
        mutual information

          Alternative metric to compare classifiers suggested by Wallach
          [1]_ and discussed by Murphy [2]_.
          
          .. note:: This metric cannot be calculated from probabilities

          Consider the following confusion matrices for classifiers A, B, and C.
          The prevalance rate is 90%.

          ::

                      A          B          C
                 ---------- ---------- ----------
                   1    0     1    0     1    0
              --+----------+----------+----------+----------+
              1 | 90   10  | 80     0 | 78     0 | HI    MI |
              0 |  0    0  | 10    10 | 12    10 | FA    CR |
              --+----------+----------+----------+----------+

          The above classifiers yield:

            ================== ======= ======= =======
            Measure               A       B       C
            ================== ======= ======= =======
            d'                  0.000   2.576   2.462
            Accuracy            0.900   0.900   0.880
            Precision           0.900   1.000   1.000
            Recall              1.000   0.888   0.867
            F-score             0.947   0.941   0.929
            Mutual information  0.000   0.187   0.174
            ================== ======= ======= =======

            Intuition suggests B > C > A but only d' and mutual information
            reflect this relationship. Mutual information is slightly more
            sensitive ([1]_ and [2]_ do discuss d').
        
          .. [1] Wallach. H. (2006). Evaluation metrics for hard classi?ers.
                 Technical report, Cavendish Lab.,Univ. Cambridge.

          .. [2] Murphy, K. P. (2007). Performance evaluation of binary
                 classifiers. http://bit.ly/LzD5m0        
        """
        # we need to determine some probabilities to aid calculating the
        # mutual information y^ refers to the predicted labels and y
        # refers to the true labels
        N = float(sum(self.values()))
        p = { 'y'  : [(self[CR]+self[FA])/N, (self[HI]+self[MI])/N],
              'y^' : [(self[CR]+self[MI])/N, (self[HI]+self[FA])/N],
             ('y^','y')  : [[self[CR]/N, self[MI]/N],
                            [self[FA]/N, self[HI]/N]] }

            
        mi = 0.
        for i,j in zip([0,0,1,1], [0,1,0,1]):
            if p[('y^','y')][i][j]:
                mi += p[('y^','y')][i][j] * \
                      math.log( p[('y^','y')][i][j] / (p['y^'][i]*p['y'][j]) )
        return mi
    
#
# Code to Implement "direct" and "prob" methods
#


# It makes sense to have a singleton so we don't have a bazzilon SDT
# instances floating around.
#
# see: http://www.garyrobinson.net/2004/03/python_singleto.html
class _S(Singleton):
    global HI,MI,CR,FA
    def __init__(self):
        self.sdt = SDT()

    def setdirects(self, hi,mi,cr,fa):
        self.sdt._directmode = True
        
        self.sdt[HI] = hi
        self.sdt[MI] = mi
        self.sdt[CR] = cr
        self.sdt[FA] = fa

    def setprobs(self, phi, pfa):
        self.sdt._directmode = False
        
        self.sdt.pHI = phi
        self.sdt.pFA = pfa
        
# let's make life a little easier by automating the
# creation of our "direct" and "prob" SDT methods


# some metrics support prob and some don't
#
# this function is dynamically loaded as a method of _vmethod
# depending on how it is initialized.
def _prob(cls, *args):
    """
    Calculates metric based on hit rate and false alarm rate
    """
    global _S
    
    if all(_isint(arg) for arg in args):
        _S.getInstance().setprobs(*args)
        func = getattr(_S.getInstance().sdt, cls.__name__)
        return func()
    else:
        results = []
        for unpacked_args in zip(*args):
            _S.getInstance().setprobs(*unpacked_args)
            func = getattr(_S.getInstance().sdt, cls.__name__)
            results.append(func())
        return results

    # Don't bothering cleaning up _S.sdt. Leave it out of
    # _directmode with whatever data happens to be there.
    # Whatever uses it next is responsible for setting it up.
    
class _vmethod(object):
    """
    Defines a factory to vectorized methods.
    """
    def __init__(self, methodname, add_prob_method=False):
        self.__name__ = methodname
        self.__doc__ = getattr(SDT, methodname).__doc__

        self._has_prob_method = add_prob_method
        if add_prob_method:
            self.prob = lambda *args: _prob(self, *args)
            self.prob.__doc__ = 'Calculates metric based on hit '\
                                'rate and false alarm rate'

    def direct(self, *args):
        """
        Calculates metric based on hit, miss, correct
        rejection, and false alarm counts
        """
        global _S
        
        if all(_isint(arg) for arg in args):
            _S.getInstance().setdirects(*args)
            func = getattr(_S.getInstance().sdt, self.__name__)
            return func()
        else:
            results = []
            for unpacked_args in zip(*args):
                _S.getInstance().setdirects(*unpacked_args)
                func = getattr(_S.getInstance().sdt, self.__name__)
                results.append(func())
            return results

        # Don't bothering cleaning up _S.sdt. Leave it in
        # _directmode with whatever data happens to be there.
        # Whatever uses it next is responsible for setting it up.

    def __call__(self, *args):
        """
        based on the number of args and the availability of .prob
        routes call to appropriate method.
        """
        if self._has_prob_method and len(args) == 2:
            return self.prob(*args)
        elif len(args) == 4:
            return self.direct(*args)
        else:
            raise Exception('Cannot route args to method')

# build methods using factory
# docstrings are grabbed from SDT methods
aprime      = _vmethod('aprime', True)
amzs        = _vmethod('amzs', True)
bpp         = _vmethod('bpp', True)
bph         = _vmethod('bph', True)
bppd        = _vmethod('bppd', True)
bmz         = _vmethod('bmz', True)
b           = _vmethod('b', True)
dprime      = _vmethod('dprime', True)
beta        = _vmethod('beta', True)
c           = _vmethod('c', True)
accuracy    = _vmethod('accuracy', True)
mcc         = _vmethod('mcc', True)
precision   = _vmethod('precision')
recall      = _vmethod('recall')
f1          = _vmethod('f1')
ppv         = _vmethod('ppv')
npv         = _vmethod('npv')
fdr         = _vmethod('fdr')
sensitivity      = _vmethod('sensitivity')
specificity      = _vmethod('specificity')
mutual_info      = _vmethod('mutual_info')
loglinear_bppd   = _vmethod('loglinear_bppd')
loglinear_dprime = _vmethod('loglinear_dprime')
loglinear_beta   = _vmethod('loglinear_beta')
loglinear_c      = _vmethod('loglinear_c')
