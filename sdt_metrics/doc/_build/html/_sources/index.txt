.. sdt_metrics documentation master file

.. currentmodule:: sdt_metrics

.. toctree::
   :hidden:

   install
   aprime
   amzs
   bpp
   bph
   bppd
   bmz
   b
   dprime
   beta
   c
   accuracy
   mcc
   precision
   recall
   f1
   ppv
   npv
   fdr
   sensitivity
   specificity
   mutual_info
   loglinear_bppd
   loglinear_dprime
   loglinear_beta
   loglinear_c
   SDT
   poc_plot
   roc_plot
   mult_roc_plot
   metric_validation_plot

:doc:`install`


Usage Examples
==============

Let's cut to the chase

    >>> from __future__ import division
    >>> from sdt_metrics import dprime
    >>> hi,mi,cr,fa = 121,42,56,34
    >>> dprime(hi,mi,cr,fa)
    0.9618717480344676

If given to arguments it'll treat them as probabilities

    >>> phi = hi/(hi+mi)
    >>> pfa = fa/(cr+fa)
    >>> dprime(phi, pfa)
    0.9618717480344676

Functions also take list-like data

    >>> dprime([.5, .6, .7],
               [.5, .5, .5])
    [0.0, 0.2533471028599986, 0.5244005132792952]

sdt_metrics.SDT is a collections.Counter-like object for storing data

    >>> from sdt_metrics import HI,MI,CR,FA, SDT
    >>> sdt_obj = SDT(HI=73,FA=32)
    >>> print(sdt_obj)
    SDT(HI=73, MI=0, CR=0, FA=32)
    >>>
    >>> sdt_obj(MI) # add a miss
    >>> print(sdt_obj)
    SDT(HI=73, MI=1, CR=0, FA=32)
    >>>
    >>> sdt_obj[CR]+=5 # add 5 correct rejections
    >>> print(sdt_obj)
    SDT(HI=73, MI=1, CR=5, FA=32)
    >>>
    >>> sdt_obj.aprime() # metrics are methods of SDT
    0.7558219178082192

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Available Metrics
=================

Parametric Metrics of Sensitivity
--------------------------------------------

.. the widths defines proportions between the columns
   added table.docutils {width:100%;} to layout.html in _templates
   to make tables fill horizontally

.. csv-table::
   :header: "",""
   :widths: 1,3
   
   :doc:`dprime`           , "d'"
   :doc:`loglinear_dprime` , "loglinear d'"
   
   
Parametric Metrics of Response Bias
--------------------------------------------

.. csv-table::
   :header: "",""
   :widths: 1,3

   :doc:`beta`             , "beta"
   :doc:`c`                , "c"
   :doc:`loglinear_beta`   , "loglinear beta"
   :doc:`loglinear_c`      , "loglinear c"
   
   
Nonparametric Metrics of Sensitivity
--------------------------------------------

.. csv-table::
   :header: "",""
   :widths: 1,3
   
   :doc:`accuracy`         , "(1 + p(HI) - p(FA)) / 2"
   :doc:`aprime`           , "A': Pollack, I., Norman, D. A. (1964)"
   :doc:`amzs`             , "Amzs: Zhang, J., and Mueller, S. T. (2005)"
   
   
Nonparametric Metrics of Response Bias
--------------------------------------------

.. csv-table::
   :header: "",""
   :widths: 1,3
   
   :doc:`b`                , "bar-napkin measure of bias: 0.5*p(HI) + 0.5p(FA)"
   :doc:`bpp`              , "B'': Grier, J. B. (1971)"
   :doc:`bph`              , "B'H: Hodos, W. (1970)"
   :doc:`bppd`             , "B''D: Donaldson, W. (1992)"
   :doc:`bmz`              , "Bmz: Zhang, J., and Mueller, S. T. (2005)"
   :doc:`loglinear_bppd`   , "loglinear B''D"
   
   
Probability Related Measures
--------------------------------------------

.. csv-table::
   :header: "",""
   :widths: 1,3
   
   :doc:`mcc`              , "Matthews Correlation Coefficient"
   :doc:`ppv`              , "Positive Predictive Value (TP / (TP + FP))"
   :doc:`npv`              , "Negative Predictive Value (TN / (TN + FN))"
   :doc:`fdr`              , "False Discovery Rate (FP / (TP + FP))"
   :doc:`sensitivity`      , "TP / (TP + FN) ``sdt.recall() <==> sdt.sensitivity()``"
   :doc:`precision`        , "TP / (TP + FP) ``sdt.precision() <==> sdt.ppv()``"
   :doc:`recall`           , "TP / (TP + FN) ``sdt.sensitivity() <==> sdt.recall()``"
   :doc:`f1`               , "Harmonic mean of precision and recall"
   :doc:`specificity`      , "TN / (TN + FP)"
   :doc:`mutual_info`      , "Wallach. H. (2006); Murphy, K. P. (2007)"

SDT Object
==================
.. csv-table::
   :header: "",""
   :widths: 1,3
   
   :doc:`SDT`      , "Data structure for holding signal detection data"


Plotting
==================
.. csv-table::
   :header: "",""
   :widths: 1,3
   
   :doc:`poc_plot`               , "Probability of Occurence Curves (POC) Plot"
   :doc:`roc_plot`               , "Receiver Operating Characteristics (ROC) Plot"
   :doc:`mult_roc_plot`          , "Multiple ROC Curves Plot"
   :doc:`metric_validation_plot` , "pcolor plot of the metric over ROC space"


Indices and tables
==================

* :ref:`genindex`
* :ref:`search`

.. * :ref:`modindex`
.. * :ref:`glossary`