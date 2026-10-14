#!/usr/bin/env python3
"""Generate plate.stl for the beamTunnel drag calibration (Zhang 2012, sec 6.4).

Plate geometry (rigid beam in side-view x-y plane):
  x in [0.010000, 0.010400]   thickness t = 4e-4 m, flow direction
  y in [0.000000, 0.008000]   beam length L = 8e-3 m, wall-normal (vertical)
  z in [-2.0e-4, +2.0e-4]     beam axis = z, width = t (matches fluid z-extent)

Origin/orientation match the FSI fluid case:
  - bottom (y=0)  = no-slip wall, beam ROOT attached
  - top    (y=H=0.01) = symmetry
  - beam root at x = L_tunnel/4 = 0.01 m  (matches paper Fig 12a)
"""

import os

t  = 4.0e-4    # plate thickness = beam b = h [m]
L  = 0.008     # beam length [m]
hz = t / 2     # half-width in z (= 2.0e-4 m)

x0, x1 = 0.0100, 0.0100 + t   # 0.010000, 0.010400
y0, y1 = 0.0,    L            # 0.000000, 0.008000
z0, z1 = -hz,   +hz           # -2.0e-4, +2.0e-4


def tri(normal, v0, v1, v2):
    n = ' '.join(str(c) for c in normal)
    verts = '\n'.join(
        f'      vertex {v[0]:.8e} {v[1]:.8e} {v[2]:.8e}' for v in (v0, v1, v2)
    )
    return f'  facet normal {n}\n    outer loop\n{verts}\n    endloop\n  endfacet\n'


faces = []

# -x face (upstream, drag face)
faces += [
    tri((-1, 0, 0), (x0, y0, z0), (x0, y0, z1), (x0, y1, z1)),
    tri((-1, 0, 0), (x0, y0, z0), (x0, y1, z1), (x0, y1, z0)),
]
# +x face (downstream)
faces += [
    tri((+1, 0, 0), (x1, y0, z0), (x1, y1, z0), (x1, y1, z1)),
    tri((+1, 0, 0), (x1, y0, z0), (x1, y1, z1), (x1, y0, z1)),
]
# -y face (root at lowerWall)
faces += [
    tri((0, -1, 0), (x0, y0, z0), (x1, y0, z0), (x0, y0, z1)),
    tri((0, -1, 0), (x1, y0, z0), (x1, y0, z1), (x0, y0, z1)),
]
# +y face (free tip)
faces += [
    tri((0, +1, 0), (x0, y1, z0), (x0, y1, z1), (x1, y1, z1)),
    tri((0, +1, 0), (x0, y1, z0), (x1, y1, z1), (x1, y1, z0)),
]
# -z and +z faces (out of 2D plane)
faces += [
    tri((0, 0, -1), (x0, y0, z0), (x0, y1, z0), (x1, y1, z0)),
    tri((0, 0, -1), (x0, y0, z0), (x1, y1, z0), (x1, y0, z0)),
]
faces += [
    tri((0, 0, +1), (x0, y0, z1), (x1, y0, z1), (x0, y1, z1)),
    tri((0, 0, +1), (x1, y0, z1), (x1, y1, z1), (x0, y1, z1)),
]

os.makedirs('constant/triSurface', exist_ok=True)
out = 'constant/triSurface/plate.stl'
with open(out, 'w') as f:
    f.write('solid plate\n')
    for face in faces:
        f.write(face)
    f.write('endsolid plate\n')

print(f'Written {out}')
print(f'  x=[{x0:.6f}, {x1:.6f}]  y=[{y0}, {y1}]  z=[{z0:.3e}, {z1:.3e}]')
print(f'  t={t:.3e} m,  L={L} m,  hz={hz:.3e} m')
