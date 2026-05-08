"""
2T3PW_force.py
==============
Python translation of 2T3PW_forceV14.nb (Wolfram Mathematica)

2 Tesla 3-Pole Wiggler (3PW) — Force calculation between upper and lower jaws.

Geometry
--------
* Period = 70 mm (wiggler period along beam direction, y-axis)
* Gap    = 32 mm (vertical aperture, z-axis)
* x-axis : horizontal transverse

Yoke structure (upper half-jaw, 1/4 model with x- and y-symmetry)
  - YokePart00  : rectangular pole + back-plate, at y = 0 … 17.5 mm
  - YokePart01  : reflection of YokePart00 about y = 17.5 mm
  - YokePart02  : copy of YokePart00 translated to y = 35 mm
  - YokeSemiCircle00 : quarter-arc connector at x=40, y=0
  - YokeSemiCircle01 : half-arc  connector at x=40, y=35

Coil structure (upper jaw)
  - Coil00 : centre racetrack coil at y=0,  j = +850 A/mm²
  - Coil01 : side   racetrack coil at y=+35, j = -260 A/mm²  (current inverted)
  - Coil02 : side   racetrack coil at y=-35, j = -260 A/mm²  (current inverted)

Lower jaw = FreeSym-copy of upper jaw rotated 180° about y-axis, currents reversed.

Expected results (from original Mathematica notebook):
  Bz(0,0,0) ≈ 2.068 T
  Fz(upper in field of lower) ≈ -4390 N  (attractive)
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import radia as rad
from math import pi
from time import time

print('RADIA version:', rad.UtiVer())


# ===========================================================================
# XC06 low-carbon steel — Froelich formula from Radia Mathematica package
# RadMatXc06[] := radMatSatIso[{1.362,0.2605,0.4917},{2118.,63.06,17.138}]
# M(H) = sum_i  ms_i * H / (ks_i + H)    [H in Oe, M in T]
# Saturation: M_sat = 1.362 + 0.2605 + 0.4917 = 2.1142 T
# ===========================================================================
_MU0 = 4.0 * pi * 1e-7   # T·m/A
_OE_TO_AM = 1000.0 / (4.0 * pi)  # 79.5775 A/m per Oe


def _M_xc06_froelich(h_oe):
    """XC06 magnetisation [T] at field h_oe [Oe] via Froelich formula."""
    return (1.362  * h_oe / (h_oe + 2118.0)
          + 0.2605 * h_oe / (h_oe + 63.06)
          + 0.4917 * h_oe / (h_oe + 17.138))


def _make_xc06():
    """Return a Radia material index for XC06 soft iron.

    Uses a tabulated M(H) curve sampled from the Froelich formula that
    matches the Mathematica RadMatXc06 definition.
    H for MatSatIsoTab is mu0*H_SI = H_Oe * 1e-4  [T].
    """
    H_oe = [0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0,
            100.0, 200.0, 500.0, 1e3, 2e3, 5e3, 1e4, 5e4, 1e5, 5e5]
    table = []
    for h in H_oe:
        h_T = h * _OE_TO_AM * _MU0   # mu0 * H_SI = H_Oe * 1e-4 T
        m_T = _M_xc06_froelich(h)
        table.append([h_T, m_T])
    return rad.MatSatIsoTab(table)


# ===========================================================================
# Model builder
# ===========================================================================
def SCW(Period, Gap):
    """Build the 2T 3-Pole Wiggler Radia model.

    Parameters
    ----------
    Period : float
        Wiggler period [mm].
    Gap    : float
        Vertical half-gap * 2 [mm] (i.e., full aperture in z).

    Returns
    -------
    dict with keys:
        'Reference' : top-level container (full magnet)
        'Yoke'      : yoke container
        'Coil'      : coil container
        'upper'     : upper jaw (yoke + coils at z > 0)
        'lower'     : lower jaw (yoke + coils at z < 0)
    """
    rad.UtiDelAll()
    mat = _make_xc06()

    # ------------------------------------------------------------------
    # Geometric parameters
    # ------------------------------------------------------------------
    YokePoleHeight        = 30.
    YokePoleWidth         = 80.
    YokePoleLength        = 0.5 * Period * (8. / 12.5)
    YokePlateHeight       = 30.
    YokePlateWidth        = 115.
    YokePlateLength       = 0.5 * Period
    YokePlateEdgeLength   = 0.35 * Period
    YokePoleEdgeLength    = YokePlateEdgeLength * (8. / 12.5)

    # Subdivisions
    YokeSubdX           = 3
    YokePoleSubdY       = 3
    YokePoleSubdZ       = 4
    YokePoleSubdRatZ    = 3.
    YokeSubdTransv      = 4
    YokeCornerSubdAz    = 5
    YokeCornerSubdRatAz = 2.5

    # Coil parameters
    CoilHeight        = YokePoleHeight - 0.5
    CoilThickness     = 0.5 * (YokePlateLength - YokePoleLength)
    # CoilEdgeThickness = 0.5 * (YokePlateEdgeLength - YokePoleEdgeLength)  # unused here
    Margin            = 0.3
    YokeCoilSpage     = 0.5
    CoilCurDens       = 850.   # A/mm²  — centre coil
    SCoilCurDens      = 260.   # A/mm²  — side coils
    CoilSubd          = 11

    # Total current [A] = current density [A/mm²] × cross-section area [mm²]
    CoilCurrentTot  = CoilCurDens  * CoilThickness * CoilHeight
    SCoilCurrentTot = SCoilCurDens * CoilThickness * CoilHeight

    # Convenience vectors
    ZeroV = [0., 0., 0.]
    Vx    = [1., 0., 0.]
    Vy    = [0., 1., 0.]
    BaseV = [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]

    # ==================================================================
    # YOKE — quarter model (x > 0, y > 0), upper half-jaw (z > 0)
    # ==================================================================
    yokeall = rad.ObjCnt([])

    # -- Pole piece: x ∈ [0, 0.5·W], y ∈ [0, 0.5·L], z ∈ [Gap/2, Gap/2+H]
    P0 = [0.25 * YokePoleWidth,
          0.25 * YokePoleLength,
          0.5  * YokePoleHeight + 0.5 * Gap]
    D0 = [0.5  * YokePoleWidth,
          0.5  * YokePoleLength,
          YokePoleHeight]
    YokePole00 = rad.ObjRecMag(P0, D0, ZeroV)
    SubdPar = [[1, 1],
               [YokePoleSubdY, 1],
               [YokePoleSubdZ, YokePoleSubdRatZ]]
    YokePole00 = rad.ObjDivMag(YokePole00, SubdPar, 'pln', BaseV,
                                'kxkykz->Numb,Frame->LabTot')

    # -- Back plate: x ∈ [0, 0.5·Pw], y ∈ [0, 0.5·Pl], z ∈ [Gap/2+H, Gap/2+H+Ph]
    P0 = [0.25 * YokePlateWidth,
          0.25 * YokePlateLength,
          0.5  * Gap + YokePoleHeight + 0.5 * YokePlateHeight]
    D0 = [0.5  * YokePlateWidth,
          0.5  * YokePlateLength,
          YokePlateHeight]
    YokePlate00 = rad.ObjRecMag(P0, D0, ZeroV)
    SubdPar = [[YokeSubdTransv, 1],
               [YokeCornerSubdAz, YokeCornerSubdRatAz],
               [1, 1]]
    # Cylindrical subdivision to follow flux lines around the yoke corner
    _near = [[0., 0.5 * YokePlateLength,
              0.5 * Gap + YokePoleHeight - 1.5],
             [1., 0., 0.]]
    _far  = [0., 0.5 * YokePlateLength,
             0.5 * Gap + YokePoleHeight + YokePlateHeight]
    _ratio = 0.5 * YokePlateLength / YokePlateHeight
    YokePlate00 = rad.ObjDivMag(YokePlate00, SubdPar, 'cyl',
                                 [_near, _far, _ratio], 'Frame->LabTot')

    # -- Semicircle connectors (iron, zero magnetisation)
    #    Polygon cross-section in the (r, z) plane relative to arc centre:
    rmin = 0.5 * YokePoleLength
    poly = [[0.001, 0.001],
            [rmin,  0.001],
            [rmin,  YokePoleHeight],
            [0.001, YokePoleHeight]]

    # Quarter-arc at x=Width/2, y=0  (phi: 0 → π/2)
    YokeSemiCircle00 = rad.ObjArcPgnMag(
        [YokePoleWidth / 2., 0.],
        'z', poly, [0., pi / 2.], 8, 'nosym', ZeroV)
    rad.ObjDivMag(YokeSemiCircle00, [1, 1, [4, 3]])

    # Half-arc  at x=Width/2, y=Period/2  (phi: -π/2 → π/2)
    YokeSemiCircle01 = rad.ObjArcPgnMag(
        [YokePoleWidth / 2., 0.5 * Period],
        'z', poly, [-pi / 2., pi / 2.], 16, 'nosym', ZeroV)
    rad.ObjDivMag(YokeSemiCircle01, [1, 1, [4, 3]])

    # Translate both arcs to z = Gap/2
    trcyla = rad.TrfTrsl([0., 0., Gap / 2.])
    rad.TrfOrnt(YokeSemiCircle00, trcyla)
    rad.TrfOrnt(YokeSemiCircle01, trcyla)

    # -- Assemble one quadrant into YokePart00, then subdivide in x
    YokePart00 = rad.ObjCnt([YokePole00, YokePlate00])
    SubdPar = [[YokeSubdX, 1], [1, 1], [1, 1]]
    YokePart00 = rad.ObjDivMag(YokePart00, SubdPar, 'pln', BaseV,
                                'Frame->LabTot')

    # -- YokePart01: reflection of YokePart00 about y = YokePlateLength / 2
    _P0 = [0., 0.5 * YokePlateLength, 0.]
    YokePart01 = rad.TrfOrnt(rad.ObjDpl(YokePart00),
                              rad.TrfPlSym(_P0, Vy))

    # -- YokePart02: translation of YokePart00 to y = Period / 2
    _P0 = [0., 0.5 * Period, 0.]
    YokePart02 = rad.TrfOrnt(rad.ObjDpl(YokePart00),
                              rad.TrfTrsl(_P0))

    # -- Collect into yokeall
    rad.ObjAddToCnt(yokeall, [YokePart00, YokePart01, YokePart02,
                               YokeSemiCircle00, YokeSemiCircle01])
    rad.MatApl(yokeall, mat)

    # -- Apply symmetry: B·n = 0 at x=0 and y=0 planes
    rad.TrfZerPerp(yokeall, ZeroV, Vx)
    rad.TrfZerPerp(yokeall, ZeroV, Vy)

    # -- Duplicate to get upper and lower jaws
    yokeupper = rad.ObjDpl(yokeall, 'FreeSym->True')
    yokelower = rad.ObjDpl(yokeall, 'FreeSym->True')
    tr = rad.TrfRot([0., 0., 0.], [0., 1., 0.], pi)   # 180° rotation about y
    rad.TrfOrnt(yokelower, tr)

    # ==================================================================
    # COILS — upper jaw
    # ==================================================================
    coilcenter = rad.ObjCnt([])
    coilsides  = rad.ObjCnt([])
    coilall    = rad.ObjCnt([])

    P00 = [0., 0., 0.5 * Gap + 0.5 * CoilHeight]
    R0  = [0.5 * YokePoleLength + YokeCoilSpage + Margin,
           0.5 * YokePoleLength - YokeCoilSpage + CoilThickness - Margin]
    D0  = [YokePoleWidth, 0.]

    # Centre coil (upper jaw, y=0)
    Coil00 = rad.ObjRaceTrk(P00, R0, D0, CoilHeight, CoilSubd,
                              CoilCurDens, 'man', 'z')

    # Side coil at y = +Period/2  (created at P00, then translated + current inverted)
    Coil01 = rad.ObjRaceTrk(P00, R0, D0, CoilHeight, CoilSubd,
                              SCoilCurDens, 'man', 'z')
    Coil01 = rad.TrfOrnt(
        rad.ObjDpl(Coil01),
        rad.TrfCmbR(rad.TrfTrsl([0., 0.5 * Period, 0.]), rad.TrfInv()))

    # Side coil at y = -Period/2  (created at P00, then translated + current inverted)
    Coil02 = rad.ObjRaceTrk(P00, R0, D0, CoilHeight, CoilSubd,
                              SCoilCurDens, 'man', 'z')
    Coil02 = rad.TrfOrnt(
        rad.ObjDpl(Coil02),
        rad.TrfCmbR(rad.TrfTrsl([0., -0.5 * Period, 0.]), rad.TrfInv()))

    # Lower-jaw coils (FreeSym copy of upper-jaw coils, rotated + current reversed)
    CCoil00 = rad.ObjDpl(Coil00, 'FreeSym->True')
    rad.TrfOrnt(CCoil00, tr)
    rad.ObjScaleCur(CCoil00, -1.)

    SCoil01 = rad.ObjDpl(Coil01, 'FreeSym->True')
    rad.TrfOrnt(SCoil01, tr)
    rad.ObjScaleCur(SCoil01, -1.)

    SCoil02 = rad.ObjDpl(Coil02, 'FreeSym->True')
    rad.TrfOrnt(SCoil02, tr)
    rad.ObjScaleCur(SCoil02, -1.)

    # Sub-containers (for display / diagnostics)
    rad.ObjAddToCnt(coilcenter, [Coil00, CCoil00])
    rad.ObjAddToCnt(coilsides,  [Coil01, SCoil01, Coil02, SCoil02])

    # Upper-jaw coil container; duplicate to make coilupper and coillower
    rad.ObjAddToCnt(coilall, [Coil00, Coil01, Coil02])
    coilupper = rad.ObjDpl(coilall, 'FreeSym->True')
    coillower = rad.ObjDpl(coilall, 'FreeSym->True')
    rad.TrfOrnt(coillower, tr)
    rad.ObjScaleCur(coillower, -1.)

    # ==================================================================
    # FINAL ASSEMBLY
    # ==================================================================
    tot      = rad.ObjCnt([])
    totupper = rad.ObjCnt([])
    totlower = rad.ObjCnt([])

    rad.ObjAddToCnt(tot,      [yokeupper, yokelower, coilupper, coillower])
    rad.ObjAddToCnt(yokeall,  [yokeupper, yokelower])
    rad.ObjAddToCnt(totupper, [yokeupper, coilupper])
    rad.ObjAddToCnt(totlower, [yokelower, coillower])

    # Update coilall to contain all coils (upper + lower, for reference)
    rad.ObjAddToCnt(coilall,
                    [Coil00, CCoil00, Coil01, SCoil01, Coil02, SCoil02])

    return {
        'Reference':       tot,
        'Yoke':            yokeall,
        'Coil':            coilall,
        'coilcenter':      coilcenter,
        'coilsides':       coilsides,
        'upper':           totupper,
        'lower':           totlower,
        'CoilCurrentTot':  CoilCurrentTot,
        'SCoilCurrentTot': SCoilCurrentTot,
    }


# ===========================================================================
# Entry point
# ===========================================================================
if __name__ == '__main__':

    pper = 70.   # period [mm]
    gg   = 32.   # gap    [mm]

    print()
    print('=' * 60)
    print(f'2T 3-Pole Wiggler  (Period={pper} mm, Gap={gg} mm)')
    print('=' * 60)

    # ---- Build model ----
    print('Building model...')
    t0 = time()
    Grp = SCW(pper, gg)
    t1 = time()
    print(f'Model built in {t1 - t0:.1f} s')

    tt = Grp['Reference']

    # ---- Solve ----
    print('Solving magnetostatic problem (tol=3e-4, max_iter=1000)...')
    t0 = time()
    re = rad.Solve(tt, 0.0003, 1000)
    t1 = time()
    print(f'Solve result : {re}')
    print(f'Solve time   : {t1 - t0:.1f} s')

    # ---- Central field ----
    bz0 = rad.Fld(tt, 'bz', [0., 0., 0.])
    print()
    print(f'Central field  Bz(0,0,0) = {bz0:.4f} T')
    print(f'  (notebook reference: ~2.068 T)')

    # ---- Self-inductance and energy (notebook: "Self Inductance and Energy") ----
    # l1 = FldEnr(coilcenter, tt) / Cur²  (effective inductance from center coils)
    # l2 = FldEnr(coilsides,  tt) / Cur²  (effective inductance from both side coil pairs)
    # Self Inductance [μH] = 2e6 * (l1 + 2*l2)
    # Coil System Energy [J] = (l1 + 2*l2) * Cur²
    Cur  = Grp['CoilCurrentTot']
    SCur = Grp['SCoilCurrentTot']
    print()
    print(f'Coil current — center: {Cur:.1f} A,  side: {SCur:.1f} A')
    print('Computing self-inductance and energy  (this may take ~1 min)...')
    t0 = time()
    rad.FldCmpPrc('PrcEnergy->1e-7')
    l1 = abs(rad.FldEnr(Grp['coilcenter'], tt, [10, 10, 10]) / Cur**2)
    l2 = abs(rad.FldEnr(Grp['coilsides'],  tt, [10, 10, 10]) / Cur**2)
    t1 = time()
    L_uH = 2.0e6 * (l1 + 2. * l2)
    E_J  = (l1 + 2. * l2) * Cur**2
    print(f'Self Inductance    : {L_uH:.4f} \u03bcH  (notebook: 0.163 \u03bcH)')
    print(f'Coil System Energy : {E_J:.1f} J  (notebook: 2033.5 J)')
    print(f'Energy time        : {t1 - t0:.1f} s')

    # ---- Force between jaws ----
    print()
    print('Computing force on upper jaw in field of lower jaw...')
    VectF = rad.FldEnrFrc(Grp['upper'], Grp['lower'], 'Fxyz', [2, 2, 2])
    print(f'Force vector [Fx, Fy, Fz] = '
          f'[{VectF[0]:.1f}, {VectF[1]:.1f}, {VectF[2]:.1f}] N')
    print(f'Fz = {VectF[2]:.1f} N')
    print(f'  (notebook reference: ~-4390 N attractive)')
