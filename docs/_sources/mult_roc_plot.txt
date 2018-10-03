mult_roc_plot
===============================================

.. currentmodule:: sdt_metrics.plotting

.. autofunction:: mult_roc_plot

Example
----------------------------------

    >>> from sdt_metrics.plotting import mult_roc_plot
    >>> sdt_obj = SDT(HI=116, MI=30, CR=323, FA=80)
    >>> sdt_probs = (.97,.22)
    >>> sdt_counts = (76,67,80,65)
    >>> mult_roc_plot((sdt_obj,  'from SDT object'),
                      (sdt_probs, 'from probs'),
                      (sdt_counts, 'from counts'),
                      isopleths='c',
                      fname='mult_roc_example.png')
    
.. image:: _static/mult_roc_example.png
    :width: 600px
    :align: center
    :height: 600px
    :alt: mult_roc_example.png
    
Example with amzs and bppd
----------------------------------

    >>> mult_roc_plot(((.91,.40), 'A'),
                      ((.76,.56), 'B'),
                      ((.84,.67), 'C'),
                      metric='amzs',
                      isopleths='bppd',
                      fname='mult_roc_example02.png')
    
.. image:: _static/mult_roc_example02.png
    :width: 600px
    :align: center
    :height: 600px
    :alt: mult_roc_example02.png