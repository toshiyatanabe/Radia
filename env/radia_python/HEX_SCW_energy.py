"""
HEX-SCW Model — Energy and Inductance
======================================
Original Mathematica notebook: HEX-SCW_model-PN_energy_TT.nb
Converted to Python.

Description
-----------
Builds a Helical-Excitation SCW (Superconducting Wiggler) magnet model in Radia
and computes the stored magnetic energy and coil self-inductance.

The model consists of:
  - 16 iron poles arranged along z with period perr = 70.14 mm
  - 17 coil-core (CC) iron pieces sitting between each pair of poles
  - 17 racetrack coils wound around the CC pieces
  - Full symmetry applied in ±z (period) and ±y (top/bottom) to model only
    the upper half and quarter-period, then mirror automatically

Key output:
  - enrj : stored energy [J], computed as abs(rad.FldEnr(ccontf, obj))
    where ccontf = coil container, obj = full magnet (coils + iron)
  - Icoil : coil current per turn [A]
  - Ltot  : total inductance [H] = 2*enrj / Icoil^2

Expected results (FEM reference):
  enrj  ≈ 64000 J   (full obj, Method 2 -- from independent FEM)
  Icoil ≈ 443 A
  Ltot  ≈ ~0.65 H  (= 2*64000/443^2)
"""

import math
import time
import radia as rad
from uti_radia_vtk import ObjDrwPyVista

print('RADIA version:', rad.UtiVer())


# ===========================================================================
# Helper: XC06 low-carbon steel material (same as Self_Inductance.py)
# ===========================================================================
def make_xc06_material():
    """Return a Radia material index for XC06 low-carbon steel."""
    mu0 = 4.0 * math.pi * 1e-7
    H   = [0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 15.0, 30.0, 60.0, 100.0, 300.0, 1000.0, 5000.0]  # Oe
    M   = [0.10, 0.25, 0.55, 0.80, 1.10, 1.35, 1.55, 1.68, 1.75, 1.78,  1.82,   1.85,   1.87]  # T
    return rad.MatSatIsoTab([[H[i] * mu0, M[i]] for i in range(len(H))])


# ===========================================================================
# Helper: rounded-rectangle cross-section polygon (2-D, in XY plane)
# ===========================================================================
def rounded_rect_profile(half_w, half_h, r_corner, n_arc=10):
    """Return list of [x, y] points tracing a rounded-rectangle cross-section.

    Parameters
    ----------
    half_w, half_h : half-dimensions of the rectangle (full sizes: 2*half_w, 2*half_h)
    r_corner       : corner radius
    n_arc          : number of arc segments per quarter circle
    Returns a closed polygon (no repeated endpoint).
    """
    pts = []
    # quarter-circle centres at (±(half_w - r), ±(half_h - r))
    for (cx_sign, cy_sign, a_start, a_end) in [
        ( 1,  1,  0,          math.pi/2),    # top-right
        (-1,  1,  math.pi/2,  math.pi),      # top-left
        (-1, -1,  math.pi,    3*math.pi/2),  # bottom-left
        ( 1, -1, -math.pi/2,  0),            # bottom-right
    ]:
        cx = cx_sign * (half_w - r_corner)
        cy = cy_sign * (half_h - r_corner)
        for k in range(n_arc):
            a = a_start + k * (a_end - a_start) / n_arc
            pts.append([cx + r_corner * math.cos(a),
                        cy + r_corner * math.sin(a)])
    return pts


# ===========================================================================
# SCW model builder
# ===========================================================================
def scw():
    """Build the HEX-SCW Radia model and return the top-level container (und)."""

    # ------ containers ------
    und      = rad.ObjCnt([])   # full undulator
    pcont    = rad.ObjCnt([])   # poles
    coilcont = rad.ObjCnt([])   # coils
    cccont   = rad.ObjCnt([])   # coil cores (CC iron)

    # ------ material ------
    matpole = make_xc06_material()

    # ------ geometry parameters [mm] ------
    nit   = 16          # number of periods (iterations)
    perr  = 70.14       # period [mm]
    perh  = perr / 2    # half-period
    gapp  = 12.0        # full gap (y) [mm]

    rcorner = 15.0      # corner radius for coil/CC cross-sections
    hcc     = 70.0 / 2  # half-height of coil core = 35 mm
    wcc     = 93.0 / 2  # half-width  of coil core = 46.5 mm
    cct     = 22.2      # coil-core thickness in z [mm]

    cout  = 10.2        # coil conductor outer radial build [mm]
    insth = 2.0         # insulation thickness [mm]

    # coil-conductor thickness along z: 14 full + 3 tapered end pieces
    coutt = [cout] * 14 + [
        (13/15) * cout,
        (7/15)  * (0.4633 / 0.35)  * cout,
        (1/15)  * (0.10666 / 0.016666) * cout,
    ]  # length = nit + 1 = 17, 1-indexed in Mathematica

    # coil axial (z) thickness: 14 full + 3 tapered
    tcoil  = 22.0       # nominal coil thickness in z [mm]
    tcoilt = [tcoil] * 14 + [tcoil, tcoil * 3/4, tcoil * 1/4]

    # current density
    jbase  = -587.0     # A/mm²
    jct    = [1.0] * 14 + [1.0, 1.0, 1.0]
    jcoil  = [jbase * j for j in jct]   # length 17, 0-indexed

    # derived
    rcorner2 = rcorner + insth + cout    # outer corner radius (pole profile)
    hpole    = hcc + insth + cout        # half-height of pole
    wpole    = wcc + insth + cout        # half-width  of pole
    polet    = 12.87                     # pole thickness in z [mm]
    wcoil    = wcc
    hcoil    = hcc

    # segmentation
    segmt  = 9    # x,y sub-divisions
    segmt2 = 2    # z sub-divisions (poles: segmt2; CC: 2*segmt2)

    # ------ pole cross-section (rounded rectangle, outer profile) ------
    # 40-segment quarter arcs matching Mathematica's Table[t, {t,0,Pi/2, 2Pi/40}]
    # (11 points per quadrant → n_arc=10 arcs per quarter)
    poletab = rounded_rect_profile(wpole, hpole, rcorner2, n_arc=10)

    # ------ CC cross-section (inner rounded rectangle) ------
    cctab = rounded_rect_profile(wcc, hcc, rcorner, n_arc=10)

    # ===========================================================
    # Build poles: i = 1..nit  (1-based, matching Mathematica)
    # ===========================================================
    for i in range(1, nit + 1):
        # ObjThckPgn(x_center, lx_thickness, [[y,z],...], axis, mag)
        # Here axis='z', so cross-section is in XY plane and extruded along Z.
        # poletab is in XY, so pass as [[x,y],…]
        p = rad.ObjThckPgn(0, polet, poletab, 'z', [0, 0, 0])

        # translate along z: center of pole
        z_pos = perh * (i - 1) + 0.5 * (polet + cct)
        rad.TrfOrnt(p, rad.TrfTrsl([0, 0, z_pos]))

        rad.ObjDrwAtr(p, [0, 1, 1], 0.001)

        # subdivide: {segmt,1} in x and y, {segmt2,1} in z, frame=LabTot
        rad.ObjDivMag(p, [[segmt, 1], [segmt, 1], [segmt2, 1]], 'Frame->LabTot')

        rad.MatApl(p, matpole)
        rad.ObjAddToCnt(pcont, [p])

    # ===========================================================
    # Build coil cores (CC iron): i = 1..nit+1
    # ===========================================================
    for i in range(1, nit + 2):
        if i <= 1:
            # half-thickness end piece at z = 0
            cc = rad.ObjThckPgn(0, cct / 2, cctab, 'z', [0, 0, 0])
            rad.TrfOrnt(cc, rad.TrfTrsl([0, 0, cct / 4]))
        else:
            # full-thickness piece
            cc = rad.ObjThckPgn(0, cct, cctab, 'z', [0, 0, 0])
            rad.TrfOrnt(cc, rad.TrfTrsl([0, 0, perh * (i - 1)]))

        rad.ObjDrwAtr(cc, [0, 1, 1], 0.001)
        rad.ObjDivMag(cc, [[segmt, 1], [segmt, 1], [2 * segmt2, 1]], 'Frame->LabTot')
        rad.MatApl(cc, matpole)
        rad.ObjAddToCnt(cccont, [cc])

    # ===========================================================
    # Build racetrack coils: i = 1..nit+1
    # ObjRaceTrk([x,y,z], [rmin,rmax], [lx,ly], h, nseg, j, 'auto', 'z')
    #   rmin = rcorner + insth
    #   rmax = rcorner + insth + coutt[i]
    #   [lx,ly] = [wcoil*2 - rcorner*2, hcoil*2 - rcorner*2]  (straight section lengths)
    #   h    = coil axial thickness (tcoilt[i] or tcoilt[i]/2 for end piece)
    #   nseg = 10
    #   j    = jcoil[i] * (-1)^i  (alternating current direction)
    # ===========================================================
    lx_coil = wcoil * 2 - rcorner * 2
    ly_coil = hcoil * 2 - rcorner * 2
    rmin_coil = rcorner + insth

    for i in range(1, nit + 2):
        rmax_i = rmin_coil + coutt[i - 1]
        j_i    = jcoil[i - 1] * ((-1) ** i)

        if i <= 1:
            # first (bottom-end) coil: half thickness, translated by tcoil/4
            h_i = tcoilt[i - 1] / 2
            coil = rad.ObjRaceTrk([0, 0, 0], [rmin_coil, rmax_i],
                                   [lx_coil, ly_coil], h_i, 10, j_i, 'auto', 'z')
            rad.TrfOrnt(coil, rad.TrfTrsl([0, 0, tcoil / 4]))

        elif i >= nit + 1:
            # last (top-end) coil: full thickness, special z position
            h_i = tcoilt[i - 1]
            coil = rad.ObjRaceTrk([0, 0, 0], [rmin_coil, rmax_i],
                                   [lx_coil, ly_coil], h_i, 10, j_i, 'auto', 'z')
            z_pos = perh * (i - 1) - (cct / 2) + (h_i / 2)
            rad.TrfOrnt(coil, rad.TrfTrsl([0, 0, z_pos]))

        else:
            # regular coil: full thickness
            h_i = tcoilt[i - 1]
            coil = rad.ObjRaceTrk([0, 0, 0], [rmin_coil, rmax_i],
                                   [lx_coil, ly_coil], h_i, 10, j_i, 'auto', 'z')
            rad.TrfOrnt(coil, rad.TrfTrsl([0, 0, perh * (i - 1)]))

        rad.ObjDrwAtr(coil, [1, 0.75, 0], 0.001)
        rad.ObjAddToCnt(coilcont, [coil])

    # ===========================================================
    # Assemble iron container (poles + CC cores) and apply symmetry
    # ===========================================================
    ironcf = rad.ObjCnt([])
    rad.ObjAddToCnt(ironcf, [cccont, pcont])

    # Translate upper half to y = hpole + gapp/2
    rad.TrfOrnt(ironcf, rad.TrfTrsl([0, hpole + gapp / 2, 0]))

    # Mirror in ±z (half-period) and ±y (top/bottom)
    rad.TrfZerPara(ironcf, [0, 0, 0], [0, 0, 1])   # z-symmetry
    rad.TrfZerPara(ironcf, [0, 0, 0], [0, 1, 0])   # y-symmetry (creates bottom half)

    # ===========================================================
    # Assemble coil container and apply symmetry
    # ===========================================================
    ccontf = rad.ObjCnt([])
    rad.ObjAddToCnt(ccontf, [coilcont])

    rad.TrfOrnt(ccontf, rad.TrfTrsl([0, hpole + gapp / 2, 0]))
    rad.TrfZerPara(ccontf, [0, 0, 0], [0, 0, 1])
    rad.TrfZerPara(ccontf, [0, 0, 0], [0, 1, 0])

    # ===========================================================
    # Top-level container
    # ===========================================================
    rad.ObjAddToCnt(und, [ccontf, ironcf])

    return und


# ===========================================================================
# Main
# ===========================================================================
print('\n' + '='*60)
print('Building HEX-SCW model...')
print('='*60)

rad.UtiDelAll()
t0 = time.time()
obj = scw()
t1 = time.time()
print(f'Model built in {t1-t0:.1f} s')
print(f'Total objects: {rad.ObjCntSize(obj)}')

# Retrieve the sub-containers (they were set as module-level names inside scw();
# here we rebuild references via the containers stored in obj).
# The returned `obj` = und = ObjCnt([ccontf, ironcf])
# We need ccontf and ironcf for the energy calculation.
# Re-build by calling scw() a second time would be expensive,
# so we pass obj as both source and destination with subdivisions.

# ------ Solve magnetization ------
print('\nSolving magnetization (tolerance=1e-6, max_iter=10000)...')
t2 = time.time()
res = rad.Solve(obj, 1e-6, 10000)
t3 = time.time()
print(f'Solve result: {res}  (time: {t3-t2:.1f} s)')

# ------ Display geometry ------
print('\nDisplaying geometry...')
ObjDrwPyVista(obj, title='HEX-SCW Model')

# ------ Energy / Inductance ------
# Because ccontf and ironcf are sub-containers of obj, we can access them
# by unpacking the top-level container.
# rad.ObjCntStuf(obj) returns [ccontf_index, ironcf_index]
contents = rad.ObjCntStuf(obj)
ccontf_ref = contents[0]   # coil container
ironcf_ref = contents[1]   # iron container

print('\n' + '='*60)
print('Computing stored energy...')
print('='*60)

# Method 1: coils in field of iron only, with subdivision [4,4,4]
print('\nMethod 1: FldEnr(ccontf, ironcf, [4,4,4])  -- iron contribution only')
t4 = time.time()
enrj_iron = abs(rad.FldEnr(ccontf_ref, ironcf_ref, [4, 4, 4]))
t5 = time.time()
print(f'  enrj (iron only) = {enrj_iron:.4f} J   (time: {t5-t4:.1f} s)')
print(f'  Expected: not yet validated against FEM')

# Method 2: coils in field of full magnet (coils + iron), no subdivision
print('\nMethod 2: FldEnr(ccontf, obj)  -- full magnet, precision-based')
t6 = time.time()
enrj_full = abs(rad.FldEnr(ccontf_ref, obj))
t7 = time.time()
print(f'  enrj (full) = {enrj_full:.4f} J   (time: {t7-t6:.1f} s)')
print(f'  Expected: ~64000 J  (independent FEM reference)'), print(f'  Note: values above ~200000 J indicate a remaining bug')

# Method 3: full magnet with subdivision [8,8,8]
#print('\nMethod 3: FldEnr(ccontf, obj, [8,8,8])  -- full magnet, subdivision')
#t8 = time.time()
#enrj_sub8 = abs(rad.FldEnr(ccontf_ref, obj, [8, 8, 8]))
#t9 = time.time()
#print(f'  enrj (seg=[8,8,8]) = {enrj_sub8:.4f} J   (time: {t9-t8:.1f} s)')
#print(f'  Expected: ~113471 J')

# ------ Inductance ------
# Icoil = J * (cross-section area) / N_turns
# From notebook: Icoil = 587 * cct * cout / 300
# Parameters: jbase=587 A/mm², cct=22.2 mm, cout=10.2 mm, 300 turns assumed
cct  = 22.2
cout = 10.2
Icoil = 587.0 * cct * cout / 300.0
print(f'\nCoil current Icoil = {Icoil:.4f} A  (expected: ~443 A)')

# Use Method 2 result for inductance (full magnet, precision-based)
enrj = enrj_full
Ltot = 2.0 * enrj / Icoil**2
print(f'Total inductance Ltot = 2*enrj/Icoil^2 = {Ltot:.6f} H')
print(f'Expected: ~0.65 H  (from FEM 64000 J reference)')
print(f'\nLtot = {Ltot*1e3:.3f} mH')
