# Copyright 2012 Patrick Varilly
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import division

import uwham
import numpy as np
import math
import scipy.stats

# Test sampling of a Gaussian...
mu_star = 4.0
sigma_star = 1.0

# ...with lots of Gaussian umbrellas
np.random.seed(1)
K = 30
sigma_k = np.random.normal(1.0, 0.25, (K,))
sigma_k[sigma_k < 0.25] = 0.25  # Ensure sigma_k isn't too small
mu_k = np.random.uniform(0.0, 10.0, (K,))
N_k = np.asarray(np.ceil(np.random.uniform(500.0, 1500.0, (K,))),
                 dtype=np.int)

# Generate samples
u_kln = np.zeros([K, K, N_k.max()], np.float64)
samples_ln = np.zeros([K, N_k.max()], np.float64)
for l in xrange(K):
    
    sample_sigma = (sigma_star ** -2 + sigma_k[l] ** -2) ** (-0.5)
    sample_mu = sample_sigma ** 2 * (mu_star / sigma_star ** 2 +
                                     mu_k[l] / sigma_k[l] ** 2)
    
    samples_ln[l, :N_k[l]] = np.random.normal(
        loc=sample_mu, scale=sample_sigma, size=[N_k[l]])
    
    for k in xrange(K):
        u_kln[k, l, :N_k[l]] = (samples_ln[l, :N_k[l]] - mu_k[k])**2 / (2*sigma_k[k]**2)

# Run UWHAM
results = uwham.UWHAM(u_kln, N_k)

# Collect statistics
dN = 0.5
x, y, y_expect = [], [], []
expect_dist = scipy.stats.norm(loc=mu_star, scale=sigma_star)
for Nmin in np.arange(0.0, 10.0, dN):
    Nmax = Nmin + dN
    P = 0.0
    
    for l in xrange(K):
        for n in xrange(N_k[l]):
            if Nmin <= samples_ln[l, n] < Nmax:
                P += results.weight_ln[l, n]

    x.append(0.5*(Nmin + Nmax))
    y.append(P)
    y_expect.append(expect_dist.cdf(Nmax) - expect_dist.cdf(Nmin))
    
x = np.asarray(x)
y = np.asarray(y)
y_expect = np.asarray(y_expect)

# Make a plot to compare
import matplotlib.pyplot as plt
plt.plot(x, np.log(y_expect), label='Exact')
plt.plot(x, np.log(y), 'ro', label='UWHAM')
plt.ylabel("P(x)")
plt.xlabel("x")
plt.legend()
plt.show()
