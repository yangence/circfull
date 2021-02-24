#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages
from circfull.version import __version__

setup(name='circfull',
    version=__version__,
    description='circfull: a tool to detect and quantify full-length circRNA isoforms from circFL-seq',
    author='Zelin Liu',
    author_email='zlliu@bjmu.edu.cn',
    url='https://github.com/yangence/circfull',
    license='GPL3',
    keywords='circular RNAs',
    python_requires=">=3",
    packages=find_packages(),
    scripts=[
        'bin/createGeneRef'
    ],
    data_files=[('bin',['bin/TRF/trf',
        'bin/TH/TideHunter'])],
    install_requires=[
        'scipy>=1.2.1',
        'pysam>=0.15.2',
        'pandas>=0.24.2',
        'numpy>=1.16.6',
        'docopt>=0.6.2',
        'python-intervals',
        'interval',
        'progressbar',
        'pyfasta',
        'mappy'
    ],
    entry_points={
      'console_scripts': [
          'circfull=circfull.circFL_main:main',
          'porechop=bin.porechop:main'
      ],
    }
)
