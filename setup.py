# Copyright (c) 2012, Roger Lew [see LICENSE.txt]
# This software is funded in part by NIH Grant P20 RR016454.

##from distutils.core import setup
from setuptools import setup


setup(name='sdt_metrics',
    version='0.1.0.0,
    description='Signal Detection Theory (SDT) metrics for Python',
    author='Roger Lew',
    author_email='rogerlew@gmail.com',
    license = "BSD",
    classifiers=["Development Status :: 5 - Production/Stable",
                 "Intended Audience :: Developers",
                 "Intended Audience :: Information Technology",
                 "Intended Audience :: Science/Research",
                 "License :: OSI Approved :: BSD License",
                 "Natural Language :: English",
                 "Programming Language :: Python :: 2.7",
                 "Topic :: Scientific/Engineering :: Bio-Informatics",
                 "Topic :: Scientific/Engineering :: Information Analysis",
                 "Topic :: Scientific/Engineering :: Mathematics",
                 "Topic :: Scientific/Engineering :: Medical Science Apps.",
                 "Topic :: Software Development :: Libraries :: Python Modules"],
    url='http://code.google.com/p/sdt-metrics/',
    packages=['sdt_metrics',
              'sdt_metrics.support',
              'sdt_metrics.analysis,
              'sdt_metrics.tests'],
    zip_safe=False)

"""C:\Python27\python.exe setup.py sdist upload --identity="Roger Lew" --sign"""
