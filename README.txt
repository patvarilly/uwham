A pure Python implementation of UWHAM.

Copyright 2012 Patrick Varilly.

Introduction
------------

UWHAM is an algorithm for stitching together data from umbrella sampling
simulations, and is closely related to MBAR and WHAM.  The key difference
is that the governing equations are solved by minimization of a convex
function instead of through self-consistent iteration, and so convergence
is *much* faster.

UWHAM:
  Zhiqiang Tan, Emilio Gallicchio, Mauro Lapelosa, and Ronald M. Levy,
  "Theory of binless multi-state free energy estimation with applications
  to protein-ligand binding", J. Chem. Phys. 136, 144102 (2012)
  http://dx.doi.org/10.1063/1.3701175 (14 pages)

MBAR:
   Michael R. Shirts and John D. Chodera,
   "Statistically optimal analysis of samples from multiple equilibrium states",
   J. Chem. Phys. 129, 124105 (2008)
   http://dx.doi.org/10.1063/1.2978177 (10 pages)

WHAM: Any modern molecular simulations textbook.  Original reference:
   Alan M. Ferrenberg and Robert H. Swendsen,
   "Optimized Monte Carlo data analysis",
   Phys. Rev. Lett. 63, 1195--1198 (1989)
   http://dx.doi.org/10.1103/PhysRevLett.63.1195 (4 pages)

Installation
------------

  The uwham package relies on NumPy and SciPy to do its work.  Make sure
  these are installed (for Mac OS X users, see note below).  Then run the
  following command inside the uwham directory:

    python setup.py install --user

Basic usage
-----------

  import uwham
  import numpy as np

  # K umbrellas, N_k independent samples in each, labeled x_kn
  # kth umbrella has biasing potential u_k(x) in units of kT
  
  N_k = [...]     # Number of samples in each umbrella
  u_kln = [...]   # u_k(x_ln), in units of kT
                  #  The value of u_kln[k, l, n] is ignored if n >= N_k[l]

  results = uwham(u_kln, N_k)

  # weight_ln[l, n] = weight of x_ln in the absence of biasing potential
  # Its value is np.NaN if n >= N_k[l]
  assert np.abs(np.nansum(results.weight_ln) - 1.0) < 1e-5

  # f_k[k] = Free energy of umbrella k with respect to unbiased ensemble,
  # in units of kT
  print results.f_k

Detailed Example
---------------

  See test_uwham.py for a self-contained example of sampling a Gaussian
  distribution with lots of parabolic umbrellas.



Installation in OS X
--------------------

In Mac OS X, the usual automatic downloading of dependencies during a
python setup.py install isn't able to successfully install NumPy and SciPy
(in my machine, it's the lack of a Fortran compiler, but the SciPy docs
point to other potential problems).  So you have to download and install
NumPy and SciPy manually.

The good news: binaries are easily available
The bad news: they only work with the version of Python from
  http://www.python.org, not the version that ships with OS X!

So you have to download and install three things:

1. Python 2.7.2 from http://www.python.org/download
2. Latest version of NumPy from http://numpy.org
3. Latest version of SciPy from http://scipy.org

This is far less painful than it sounds.  In my own case, the disk images
that I ended up downloading were:

1. python-2.7.2-macosx10.6.dmg
2. numpy-1.6.1-py2.7-python.org-macosx10.6.dmg
3. scipy-0.10.0-py2.7-python.org-macosx10.6.dmg

CAREFUL: for some of these packages, the "link to the latest version" that
SourceForge suggests may be incorrect!  Do look at the full list of downloads
available and pick the one that is most appropriate to your own setup

Once you've installed python2.7, numpy and scipy, you can run the command

    python setup.py install --user
