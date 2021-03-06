import sdt_metrics
from string import Template

template = Template("""\
.. This file was autogenerated by the build_metric_rst_files.py script

$metric
===============================================

.. currentmodule:: sdt_metrics


Download 
[:download:`hi-res <_static/${metric}.png>`] 
[:download:`pdf <_static/${metric}.pdf>`] 

.. image:: _static/${metric}__lores.png 
    :width: 700px
    :align: center
    :height: 341px
    :alt: ${metric}__lores.png
    
.. autoclass:: sdt_metrics.$metric
   :members: $members""")

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
    print(metric)

    if hasattr(getattr(sdt_metrics, metric), 'prob'):
        members = 'direct, prob, __call__'
    else:
        members = 'direct, __call__'

##    if metric in ['beta','bmz']:
##        validation_plot = 'log(%s).png'%metric
##    elif metric is 'bph':
##        validation_plot = 'sgn(bph)log(abs(bph+1)).png'
##    else:
##        validation_plot = '%s.png'%metric
                
    with open('%s.rst'%metric, 'wb') as f:
        f.write(template.substitute(metric=metric,
                                    members=members))
    
