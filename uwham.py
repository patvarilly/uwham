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

"""A pure Python implementation of UWHAM.

Reference:

 Zhiqiang Tan, Emilio Gallicchio, Mauro Lapelosa, and Ronald M. Levy,
 "Theory of binless multi-state free energy estimation with applications
 to protein-ligand binding", J. Chem. Phys. 136, 144102 (2012)
 http://dx.doi.org/10.1063/1.3701175 (14 pages)
"""

__all__ = ['UWHAM']

import numpy as np
import scipy.optimize

class UWHAM(object):
    """Reweight all samples from a set of umbrella sampling simulations.

    Suppose you perform K umbrella sampling simulations.  The kth simulation
    (k = 0, ..., K-1) is done with a reduced biasing potential u_k(x) (i.e.,
    u_k(x) is in units of kT), and you obtain N_k independent samples,
    labelled x_kn (n = 0, ..., N_k - 1).  As long as there is sufficient
    overlap between the umbrellas, this code will compute the weight
    that each sample x_kn would have in the absence of any biasing
    potential.  Conceptually, it is similar to WHAM and MBAR, but the
    equations that govern the problem are solved by minimizing a convex
    function instead of through self-consistent iteration.  Hence,
    convergence is *much* faster.

    Parameters
    ----------
    u_kln : 3D array, dimensions K x K x max(N_k)
      u_kln[k, l, n] is equal to u_k(x_ln).  If n exceed N_l, the value
      of u_kln is ignored.

    N_k : 1D array, dimension K
      N_k[k] is the number of independent samples in umbrella k.

    Attributes
    ----------
    weight_ln : 2D array, dimensions K x max(N_k)
      weight_ln[l, n] is the weight of the sample x_ln in the absence of
      any biasing potential.  If n exceeds N_l, the weight is NaN.
      Note that np.nansum(self.weight_ln) == 1.0.
      
    f_k : 1D array, dimension K
      f_k[k] = free energy difference (in units of kT) between umbrella k
      and the unbiased ensemble (i.e. f_k = -ln(Z_k / Z_0)).
    """

    def __init__(self, u_kln, N_k):
        # Make copies of input arrays (so we can modify u_kln below)
        self.N_k = np.array(N_k)
        self.u_kln = np.array(u_kln)

        self.K = N_k.shape[0]

        self.f_k = np.zeros([self.K], np.float64)
        self.ln_N_k = np.log(N_k)

        # Set all nonsense u_kln to NaN
        for k, l in np.ndindex(self.u_kln.shape[0:2]):
            self.u_kln[k, l, N_k[l]:] = np.NaN

        # Do it
        self._solve()

    def _ln_sum_exp(self, A, axis=0):
        """Calculate ln( sum_i exp( A_i ) ), where the index runs
        along the given axis of A."""
        max_A = np.nanmax(A, axis=axis)
        return np.log(np.nansum(np.exp(
                    A - np.expand_dims(max_A, axis)), axis=axis)) + max_A
        
    def _recalc_lnW_ln(self):
        # lnW_ln = ln[ sum_k N_k exp{f_k - u_k(x_ln)} ]
        self._lnW_ln = self._ln_sum_exp(
            (self.ln_N_k + self.f_k)[:, np.newaxis, np.newaxis] - self.u_kln,
            axis=0)

    def _solve(self):
        # Function kappa to minimize, with the constraint that f_k[0] == 0.0
        def kappa(other_f_ks):
            
            if any(self.f_k[1:] != other_f_ks):
                self.f_k[0] = 0.0
                self.f_k[1:] = other_f_ks
                self._recalc_lnW_ln()
            
            return np.nansum(self._lnW_ln) - np.sum(self.N_k * self.f_k)

        def grad_kappa(other_f_ks):
            
            if any(self.f_k[1:] != other_f_ks):
                self.f_k[0] = 0.0
                self.f_k[1:] = other_f_ks
                self._recalc_lnW_ln()
                
            return (
                np.exp(
                    self._ln_sum_exp(
                        np.reshape(
                            (self.ln_N_k + self.f_k)[:, np.newaxis, np.newaxis]
                            - self.u_kln
                            - self._lnW_ln[np.newaxis, :, :],
                            (self.K, -1)),
                        axis=1))
                - self.N_k)[1:]

        def callback(other_f_ks):
            print other_f_ks

        # Minimize kappa
        self.f_k[:] = np.NaN
        
        self.f_k[1:] = scipy.optimize.fmin_bfgs(
            f=kappa,
            x0=np.zeros([self.K - 1], np.float64),
            fprime=grad_kappa,
            maxiter=1000,
            callback=callback)
        
        self.f_k[0] = 0.0
        self._recalc_lnW_ln()

        # Now renormalize free energies to have exp(-lnW_ln) be the
        # normalized probability of sample x_ln in the absence of a bias
        f_unbiased = -self._ln_sum_exp(-self._lnW_ln.ravel())
        self.f_k -= f_unbiased
        self._lnW_ln -= f_unbiased

        self.weight_ln = np.exp(-self._lnW_ln)
