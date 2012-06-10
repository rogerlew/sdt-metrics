roc_plot
===============================================

.. currentmodule:: sdt_metrics.plotting

.. autofunction:: roc_plot

Producing plot from probabilities
----------------------------------

    >>> from sdt_metrics.plotting import roc_plot
    >>> roc_plot(.67, .43, metric='amzs', fname='roc_example01.png')
    
.. image:: _static/roc_example01.png
    :width: 400px
    :align: center
    :height: 400px
    :alt: roc_example01.png
    
Producing plot from frequencies
----------------------------------

    >>> from sdt_metrics.plotting import roc_plot
    >>> roc_plot(134, 23, 345, 67, fname='roc_example01.png')
    
.. image:: _static/roc_example02.png
    :width: 400px
    :align: center
    :height: 400px
    :alt: roc_example02.png
    
Producing plot from :class:`SDT` Object
----------------------------------------

    >>> from sdt_metrics import HI,MI,CR,FA, SDT
    >>> from random import choice
    >>> sdt_obj = SDT([choice([HI,MI,CR,FA]) for i in xrange(1000)])
    >>> print(sdt_obj)
    SDT(HI=251, MI=245, CR=264, FA=240)
    >>> roc_plot(sdt_obj, fname='roc_example03.png')
    
.. image:: _static/roc_example03.png
    :width: 400px
    :align: center
    :height: 400px
    :alt: roc_example03.png