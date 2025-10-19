#!/usr/bin/env python3
"""Description:

Setup script for WANG

"""
import platform

""" The following codes were copied from MACS project and will be modified in the future."""
import sys
import os
import re
from setuptools import setup, Extension
from Cython.Build import cythonize
import subprocess
import sysconfig
import numpy

# get version
exec(open("WangLab/Constants.py").read())

# classifiers
classifiers = [ \
    'Development Status :: 5 - Production/Stable',
    'Environment :: Terminal',
    'Intended Audience :: Developers',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: MIT License',
    'Operating System :: MacOS :: MacOS X',
    'Operating System :: POSIX',
    'Operating System :: Unix',
    'Operating System :: Windows',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'Programming Language :: Python :: 3.9',
    'Programming Language :: Python :: 3.10',
    'Programming Language :: Python :: 3.11',
    'Programming Language :: Python :: 3.12',
    'Programming Language :: Cython', ]

install_requires = ["argcomplete",
                    "biopython",
                    "pandas",
                    "numpy",
                    "selenium",
                    "cutadapt",
                    "openpyxl",
                    ]
tests_requires = ['pytest']
scripts = ['bin/wanglab', 'bin/wanglab.bat']
def main():
    if sys.version_info < (3, 9):
        sys.stderr.write("CRITICAL: Python version must >= 3.9!\n")
        sys.exit(1)

    # NumPy include dir
    numpy_include_dir = [numpy.get_include()]

    with open("README.md", "r") as fh:
        long_description = fh.read()

    setup(name="WangLab",
          version=VERSION,
          description="Simple bioinformatic tools",
          long_description=long_description,
          long_description_content_type="text/markdown",
          author='Yifan Bu',
          author_email='y30210580@163.com',
          url='http://github.com/byf1999/WangLab/',
          package_dir={'WangLab': 'WangLab'},
          packages=['WangLab', 'WangLab.ChIP_seq', 'WangLab.RNA_seq', 'WangLab.Sequence_operate', 'WangLab.TIS'],
          package_data={'WangLab': ['*.pxd']},
          scripts=scripts,
          classifiers=classifiers,
          install_requires=install_requires,
          tests_require=tests_requires,
          python_requires='>=3.9',
          )


if __name__ == '__main__':
    main()
