##############################################################################
# RADIA Insertion Device (ID) Utility Functions
# Provides solve and field-calculation helpers used by SCU.py and similar scripts.
##############################################################################

from __future__ import absolute_import, division, print_function
import radia as rad
import sys
import time
from math import sqrt

# ---------------------------------------------------------------------------
# Solve (relaxation) helper
# ---------------------------------------------------------------------------

def id_solve(_id, _prc=0.001, _max_it=5000, _print=2):
    """Build the interaction matrix and run the relaxation for a Radia object."""
    if _print > 0:
        print('Setting up interaction matrix ...')
        sys.stdout.flush()
    t0 = time.time()
    IM = rad.RlxPre(_id)
    if _print > 0:
        print('Interaction matrix set up in: {:.2f} s'.format(time.time() - t0))
        sys.stdout.flush()

    if _print > 0:
        print('Performing relaxation (tol={}, max_it={}) ...'.format(_prc, _max_it))
        sys.stdout.flush()
    t0 = time.time()
    res = rad.RlxAuto(IM, _prc, _max_it)
    if _print > 0:
        print('Relaxation took: {:.2f} s'.format(time.time() - t0))
        print('Relaxation result:', res)
        sys.stdout.flush()
    return res


# ---------------------------------------------------------------------------
# Magnetic field vs. longitudinal position
# ---------------------------------------------------------------------------

def mag_fld_vs_long_pos(_id, _ns_per=50, _per=14., _nper=4, _long_ax='z',
                        _rank=0, _print=2, _half_nper_ext=8,
                        _e_beam_GeV=7.0):
    """Calculate horizontal and vertical on-axis magnetic field along an undulator.

    Returns
    -------
    arBh   : list of horizontal field values [T]
    arBv   : list of vertical field values [T]
    ns     : total number of calculation points
    sRange : longitudinal calculation range [mm]
    Keff   : effective undulator deflection parameter (dimensionless)
    E1     : estimated first-harmonic photon energy [eV]
    """
    nTotPer = _nper + 2 * _half_nper_ext
    ns = _ns_per * nTotPer + 1
    sRange = nTotPer * _per          # mm
    sStep  = sRange / (ns - 1)
    sStart = -0.5 * sRange

    if _long_ax == 'y':
        pts = [[0, sStart + i * sStep, 0] for i in range(ns)]
        Bh_comp = 'bx'
        Bv_comp = 'bz'
    else:  # default: longitudinal axis = z
        pts = [[0, 0, sStart + i * sStep] for i in range(ns)]
        Bh_comp = 'bx'
        Bv_comp = 'by'

    t0 = time.time()
    arBh = rad.Fld(_id, Bh_comp, pts)
    arBv = rad.Fld(_id, Bv_comp, pts)

    if _print > 0 and _rank <= 0:
        print('Field calculation took: {:.2f} s'.format(time.time() - t0))
        sys.stdout.flush()

    # Peak fields
    if isinstance(arBh, (list, tuple)):
        B0h = max(abs(b) for b in arBh)
        B0v = max(abs(b) for b in arBv)
    else:
        B0h = abs(arBh)
        B0v = abs(arBv)

    # K parameter:  K = 0.09336 * B0[T] * lambda_u[mm]
    Kh   = 0.09336 * B0h * _per
    Kv   = 0.09336 * B0v * _per
    Keff = sqrt(Kh**2 + Kv**2)

    # First-harmonic photon energy [eV]:
    #   E1 = 0.9496e6 * E_GeV^2 / ((1 + K^2/2) * lambda_u[cm])
    lam_u_cm = _per * 0.1
    E1 = 0.9496e6 * _e_beam_GeV**2 / ((1. + Keff**2 / 2.) * lam_u_cm)

    if _print > 0 and _rank <= 0:
        print('Peak vertical field: B0v = {:.4f} T'.format(B0v))
        print('Keff = {:.4f}'.format(Keff))
        print('E1   = {:.1f} eV'.format(E1))
        sys.stdout.flush()

    return arBh, arBv, ns, sRange, Keff, E1
