#!/usr/bin/env python3
"""Generate plate.stl for sideView drag calibration.

Plate geometry (full beam cross-section in x-y plane):
  x ∈ [0.020000, 0.020212]  (thickness t = 2.12e-4 m, flow direction)
  y ∈ [0.000000, 0.008000]  (beam length L = 8e-3 m, wall-normal)
  z ∈ [-1.06e-4, +1.06e-4]  (beam axis = z, width = t)
"""

import os

t  = 2.12e-4   # plate thickness = beam b = h [m]
L  = 0.008     # beam length [m]
hz = t / 2     # half-width in z (= 1.06e-4 m)

x0, x1 = 0.02,  0.02 + t   # 0.020000, 0.020212
y0, y1 = 0.0,   L          # 0.000000, 0.008000
z0, z1 = -hz,  +hz         # -1.06e-4, +1.06e-4

def tri(normal, v0, v1, v2):
    n = ' '.join(str(c) for c in normal)
    verts = '\n'.join(
        f'      vertex {v[0]:.8e} {v[1]:.8e} {v[2]:.8e}' for v in (v0, v1, v2)
    )
    return f'  facet normal {n}\n    outer loop\n{verts}\n    endloop\n  endfacet\n'

faces = []

# -x face (x=x0, outward normal (-1,0,0)): upstream drag face
faces += [
    tri((-1,0,0), (x0,y0,z0), (x0,y0,z1), (x0,y1,z1)),
    tri((-1,0,0), (x0,y0,z0), (x0,y1,z1), (x0,y1,z0)),
]
# +x face (x=x1, outward normal (+1,0,0)): downstream face
faces += [
    tri((+1,0,0), (x1,y0,z0), (x1,y1,z0), (x1,y1,z1)),
    tri((+1,0,0), (x1,y0,z0), (x1,y1,z1), (x1,y0,z1)),
]
# -y face (y=y0=0, outward normal (0,-1,0)): plate base at lowerWall
faces += [
    tri((0,-1,0), (x0,y0,z0), (x1,y0,z0), (x0,y0,z1)),
    tri((0,-1,0), (x1,y0,z0), (x1,y0,z1), (x0,y0,z1)),
]
# +y face (y=y1=L, outward normal (0,+1,0)): plate free tip
faces += [
    tri((0,+1,0), (x0,y1,z0), (x0,y1,z1), (x1,y1,z1)),
    tri((0,+1,0), (x0,y1,z0), (x1,y1,z1), (x1,y1,z0)),
]
# -z face (z=z0, outward normal (0,0,-1))
faces += [
    tri((0,0,-1), (x0,y0,z0), (x0,y1,z0), (x1,y1,z0)),
    tri((0,0,-1), (x0,y0,z0), (x1,y1,z0), (x1,y0,z0)),
]
# +z face (z=z1, outward normal (0,0,+1))
faces += [
    tri((0,0,+1), (x0,y0,z1), (x1,y0,z1), (x0,y1,z1)),
    tri((0,0,+1), (x1,y0,z1), (x1,y1,z1), (x0,y1,z1)),
]

os.makedirs('constant/triSurface', exist_ok=True)
out = 'constant/triSurface/plate.stl'
with open(out, 'w') as f:
    f.write('solid plate\n')
    for face in faces:
        f.write(face)
    f.write('endsolid plate\n')

print(f'Written {out}')
print(f'  x=[{x0:.6f}, {x1:.8e}]  y=[{y0}, {y1}]  z=[{z0:.3e}, {z1:.3e}]')
print(f'  t={t:.3e} m,  L={L} m,  hz={hz:.3e} m')
