#!/usr/bin/env python
#
# Setup script for uwham module.
# Written by Patrick Varilly, 5-6 May 2012

import sys

if sys.version_info[0] >= 3:
    print "WARNING: uwham package has not been tested on Python 3"

CLASSIFIERS = """\
Development Status :: 3 - Alpha
Intended Audience :: Science/Research
License :: OSI Approved :: GNU General Public License v3 (GPLv3)
Operating System :: POSIX
Operating System :: Unix
Operating System :: MacOS
Operating System :: Microsoft :: Windows
Programming Language :: Python
Topic :: Scientific/Engineering

"""

try:
    from setuptools import setup, Extension
except ImportError:
    from distutils.core import setup, Extension

# On OS X, the automatic downloading and installation of numpy and scipy
# is very problematic.  Avoid it
if sys.platform == 'darwin':
    try:
        import numpy
        import scipy
    except ImportError:
        print """

FATAL ERROR: NumPy and SciPy both need to be installed to install uwham.
Ordinarily, they would be installed automatically, but this seems to be
very problematic in OS X.  Fortunately, there are easy-to-install binaries
available for download.

Please read the README.txt file for installation instructions.

"""
        exit()

setup(name="uwham",
      version="1.0",
      author="Patrick Varilly",
      author_email="pv271@cam.ac.uk",
      description="Pure Python implementation of UWHAM (cf. MBAR, WHAM)",
      url="http://github.com/patvarilly/uwham",
      download_url="http://github.com/patvarilly/uwham/downloads",
      license='GPLv3',
      classifiers=[_ for _ in CLASSIFIERS.split('\n') if _],
      platforms=["Linux", "Unix", "Mac OS-X", "Windows", "Solaris"],
      install_requires=[
        'numpy>=1.3.0',
        'scipy>=0.7.0',
        ],
      
      py_modules=['uwham'],
      )
