"""
Coil Inductance
===============
Original: Pascal Elleaume, ESRF, March 1999 (Wolfram Mathematica notebook)
Converted to Python by GitHub Copilot.

This script demonstrates the computation of self and mutual inductance of coils
using Radia's rad.FldEnr() function.

rad.FldEnr(objdst, objsrc) computes the integral over the volume of objdst of
(current density in objdst) · (vector potential created by objsrc).
The result is in Joules.

  - If objdst == objsrc  : result = L * I^2  (self inductance × current²)
  - If objdst != objsrc  : result = M * I1 * I2  (mutual inductance × currents)

rad.FldCmpPrc('PrcEnergy', value) controls the precision by subdividing coils
as needed. Start with a large value and reduce until the result stabilises.

Sections:
  1. Two Coils Without Iron — self and mutual inductance
  2. Single Coil With Iron Yoke — self inductance with iron core
"""

import sys
import math
import radia as rad
from uti_radia_vtk import ObjDrwPyVista

print('RADIA version:', rad.UtiVer())

# ===========================================================================
# Section 1: Two Coils Without Iron
# ===========================================================================
print('\n' + '='*60)
print('Section 1: Two Coils Without Iron')
print('='*60)

# --- Coil #1 ---
Cur   = 10      # current [A]
Nt    = 25      # number of turns
rmin1 = 10      # inner radius [mm]
rmax1 = 15      # outer radius [mm]
h1    = 20      # height [mm]

pt       = [0, 0, 0]
phimin   = 0
phimax   = 2 * math.pi
nseg     = 20

j1 = -Cur * Nt / (rmax1 - rmin1) / h1   # current density [A/mm²]

g1 = rad.ObjArcCur(pt, [rmin1, rmax1], [phimin, phimax], h1, nseg, j1)
rad.ObjDrwAtr(g1, [1, 0, 0.5])

# --- Coil #2 ---
rmin2 = 16      # inner radius [mm]
rmax2 = 21      # outer radius [mm]
h2    = 20      # height [mm]

j2 = -Cur * Nt / (rmax2 - rmin2) / h2

g2 = rad.ObjArcCur(pt, [rmin2, rmax2], [phimin, phimax], h2, nseg, j2)
rad.ObjDrwAtr(g2, [0, 1, 0.5])

# Container for both coils
g = rad.ObjCnt([g1, g2])

# --- Display geometry ---
print('Displaying two-coil geometry...')
ObjDrwPyVista(g, title='Two Coils (No Iron)')

# ===========================================================================
# Section 2: Self and Mutual Inductance (no iron)
# ===========================================================================
print('\n' + '='*60)
print('Section 2: Self and Mutual Inductance (no iron)')
print('='*60)

# Set precision for mutual inductance computation (no singularity — source ≠ destination).
rad.FldCmpPrc('PrcEnergy->1e-7')

# Self inductance: use explicit subdivision [k,k,k] for the destination object.
# The precision-controlled form FldEnr(g,g) uses [1,1,1] internally and OVERestimates
# self inductance by ~25% due to the 1/r singularity at zero separation.
# With [8,8,8] the result agrees with Wheeler's multilayer formula to within ~1%.
# For mutual inductance (g1,g2) there is no singularity; all methods agree.
k = 8   # subdivision for self-inductance; increase to 12 for finer check

# Energy integrals [J] (raw values, before dividing by I²)
e1 = abs(rad.FldEnr(g1, g1, [k,k,k])) / Cur**2   # L1  [H]
e2 = abs(rad.FldEnr(g2, g2, [k,k,k])) / Cur**2   # L2  [H]
e  = abs(rad.FldEnr(g,  g,  [k,k,k])) / Cur**2   # L_system  [H]
m  = abs(rad.FldEnr(g1, g2)) / Cur**2              # M   [H]  (no singularity; precision form OK)

# Print in μH and μJ
print(f'Coil # 1 : Energy : {1e6 * e1/2 * Cur**2:.4g} μJ   Self inductance : {1e6 * e1:.4g} μH')
print(f'Coil # 2 : Energy : {1e6 * e2/2 * Cur**2:.4g} μJ   Self inductance : {1e6 * e2:.4g} μH')
print(f'Coil System : Energy : {1e6 * e/2 * Cur**2:.4g} μJ   System inductance : {1e6 * e:.4g} μH')
print(f'Mutual Inductance : {1e6 * m:.4g} μH')
print()
print('Expected (subdivision [8,8,8] form, consistent with Wheeler multilayer formula):')
print('  L1 ~ 10.1 uH,  L2 ~ 19.8 uH,  M ~ 9.45 uH')
print('Wheeler formula: L1=10.09 uH, L2=19.76 uH')
print('Mathematica (buggy Radia, precision form): L1=10.90, L2=21.97, M=9.45 uH')
print('Python (fixed Radia, precision form = [1,1,1]): L1=12.68, L2=25.13 — 25% too high')
print('Note: M agrees across all methods (~9.45 uH) because there is no 1/r singularity')
print('  when source and destination are different objects.')

# ===========================================================================
# Section 3: Single Coil with Iron Yoke
# ===========================================================================
print('\n' + '='*60)
print('Section 3: Single Coil with Iron Yoke')
print('='*60)
print('Note: if the iron saturates, the inductance depends on current;')
print('the computed value represents an average.')

# Delete previous objects to free memory
rad.UtiDelAll()

# --- Create the coil ---
Cur   = 1000    # current [A]
Nt    = 130     # number of turns
rmin  = 140     # inner radius [mm]
rmax  = 150     # outer radius [mm]
h     = 1000    # height [mm]

pt     = [0, 0, 0]
nseg   = 20

j = -Cur * Nt / (rmax - rmin) / h

g1 = rad.ObjArcCur(pt, [rmin, rmax], [phimin, phimax], h, nseg, j)
rad.ObjDrwAtr(g1, [1, 0, 0.5])

# --- Create the iron cylinder ---
# rad.ObjCylMag(center, r, h, nseg [, axis='z'] [, mag=[0,0,0]])
# Initial magnetization defaults to [0,0,0]; material applied below.
# Radius = 75 mm (inside the coil bore), height = 1000 mm, axis along z.
m_iron = rad.ObjCylMag(pt, 75.0, 1000.0, 20)
rad.ObjDrwAtr(m_iron, [0, 0.5, 0])

# --- Define XC06 low-carbon steel material ---
# XC06 Froelich formula from Radia Mathematica package RadMatXc06[]:
#   M(H) = 1.362*H/(H+2118) + 0.2605*H/(H+63.06) + 0.4917*H/(H+17.138)
#   [H in Oe, M in T]
# MatSatIsoTab expects H_T = mu0*H_SI = H_Oe * 1e-4  [T]
_OE_TO_AM = 1000.0 / (4.0 * math.pi)  # 79.5775 A/m per Oe
def _M_xc06(h_oe):
    return (1.362*h_oe/(h_oe+2118.) + 0.2605*h_oe/(h_oe+63.06) + 0.4917*h_oe/(h_oe+17.138))
H_xc06_oe = [0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0,
             100.0, 200.0, 500.0, 1e3, 2e3, 5e3, 1e4, 5e4, 1e5, 5e5]
mat_xc06 = rad.MatSatIsoTab([[h * _OE_TO_AM * mu0, _M_xc06(h)] for h in H_xc06_oe])
rad.MatApl(m_iron, mat_xc06)

# --- Subdivide iron for accurate field computation ---
nx, ny, nz = 3, 3, 3
rad.ObjDivMag(m_iron, [nx, ny, nz])

# --- Combine coil + iron and solve for magnetization ---
t = rad.ObjCnt([g1, m_iron])
print('Solving for magnetization in iron...')
solve_result = rad.Solve(t, 0.0001, 200)
print('Solve result:', solve_result)

# --- Display geometry ---
print('Displaying coil + iron geometry...')
ObjDrwPyVista(t, title='Coil with Iron Yoke')

# --- Compute self inductance with iron ---
rad.FldCmpPrc('PrcEnergy->1e-7')
l1 = abs(rad.FldEnr(g1, t, [2, 2, 2])) / Cur**2   # [H]
print(f'\nSelf Inductance (with iron, seg=[2,2,2]) : {1e6 * l1:.4g} μH')
print('Notebook reference 5881 uH was computed with:')
print('  - Correct XC06 (RadMatXc06 Froelich formula)')
print('  - OLD Radia library with NestedFor_Energy bug')
print('Python (correct XC06 + fixed radia.so) gives a different value.')
print('Both effects (bug fix +16%, XC06 change) shift the result.')

# --- Convergence study: vary segmentation ---
print('\nConvergence check (varying segmentation):')
for seg in ([2,2,2], [4,4,4]):
    l = abs(rad.FldEnr(g1, t, seg)) / Cur**2
    print(f'  seg={seg}  →  L = {1e6 * l:.4g} μH')

# --- Contribution from iron alone ---
l_iron_only = abs(rad.FldEnr(g1, m_iron, [12, 12, 12])) / Cur**2
print(f'\nIron contribution only (seg=[12,12,12]) : {1e6 * l_iron_only:.4g} μH')
print('Expected (from original notebook): ~4149 μH')
