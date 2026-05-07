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

# Set precision for energy / inductance computation.
# Reduce PrcEnergy until the result is stable; smaller values = better accuracy
# but higher memory use. Start large and work down.
rad.FldCmpPrc('PrcEnergy->1e-7')

# Energy integrals [J] (raw values, before dividing by I²)
# Dividing rad.FldEnr(a, b) by I² gives L (if a==b) or M (if a!=b) in [H]
e1 = abs(rad.FldEnr(g1, g1)) / Cur**2   # L1  [H]
e2 = abs(rad.FldEnr(g2, g2)) / Cur**2   # L2  [H]
e  = abs(rad.FldEnr(g,  g))  / Cur**2   # L_system  [H]
m  = abs(rad.FldEnr(g1, g2)) / Cur**2   # M   [H]

# Print in μH and μJ
print(f'Coil # 1 : Energy : {1e6 * e1/2 * Cur**2:.4g} μJ   Self inductance : {1e6 * e1:.4g} μH')
print(f'Coil # 2 : Energy : {1e6 * e2/2 * Cur**2:.4g} μJ   Self inductance : {1e6 * e2:.4g} μH')
print(f'Coil System : Energy : {1e6 * e/2 * Cur**2:.4g} μJ   System inductance : {1e6 * e:.4g} μH')
print(f'Mutual Inductance : {1e6 * m:.4g} μH')
print()
print('Expected (from original notebook):')
print('  Coil #1 self  : ~10.9 μH')
print('  Mutual        : ~9.45 μH')

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
# M-H curve for XC06 steel (H in Oersted, M in Tesla).
# convH = μ₀ converts H[Oe] → [T] for MatSatIsoTab input: [[μ₀·H, M], ...]
# Replace with your own measured XC06 data for higher accuracy.
mu0 = 4.0 * math.pi * 1e-7   # [T·m/A]
H_xc06 = [0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 15.0, 30.0, 60.0, 100.0, 300.0, 1000.0, 5000.0]  # Oe
M_xc06 = [0.10, 0.25, 0.55, 0.80, 1.10, 1.35, 1.55, 1.68, 1.75, 1.78,  1.82,   1.85,   1.87]  # T
mat_xc06 = rad.MatSatIsoTab([[H_xc06[i] * mu0, M_xc06[i]] for i in range(len(H_xc06))])
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
print('Expected (from original notebook): ~5881 μH')

# --- Convergence study: vary segmentation ---
print('\nConvergence check (varying segmentation):')
for seg in ([2,2,2], [4,4,4]):
    l = abs(rad.FldEnr(g1, t, seg)) / Cur**2
    print(f'  seg={seg}  →  L = {1e6 * l:.4g} μH')

# --- Contribution from iron alone ---
l_iron_only = abs(rad.FldEnr(g1, m_iron, [12, 12, 12])) / Cur**2
print(f'\nIron contribution only (seg=[12,12,12]) : {1e6 * l_iron_only:.4g} μH')
print('Expected (from original notebook): ~4149 μH')
