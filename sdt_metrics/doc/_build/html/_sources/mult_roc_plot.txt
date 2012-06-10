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
                      fname = 'mult_roc_example.png')
    
.. image:: _static/mult_roc_example.png
    :width: 400px
    :align: center
    :height: 400px
    :alt: mult_roc_example.png