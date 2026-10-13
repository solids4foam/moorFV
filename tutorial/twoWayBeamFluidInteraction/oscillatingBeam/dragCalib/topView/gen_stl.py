#!/usr/bin/env python3
"""Generate plate.stl: a box representing the square beam cross-section."""

import os

t = 2.12e-4  # plate thickness (m)

# Box corners: x in [-t/2, t/2], z in [-t/2, t/2], y in [0, t]
x0, x1 = -t/2, t/2
y0, y1 = 0.0,   t
z0, z1 = -t/2, t/2

def tri(normal, v0, v1, v2):
    n = ' '.join(f'{c}' for c in normal)
    verts = '\n'.join(
        f'      vertex {v[0]:.8e} {v[1]:.8e} {v[2]:.8e}' for v in (v0, v1, v2)
    )
    return f'  facet normal {n}\n    outer loop\n{verts}\n    endloop\n  endfacet\n'

faces = []
# -x face, outward normal (-1, 0, 0): go +z then +y
faces += [
    tri((-1,0,0), (x0,y0,z0), (x0,y0,z1), (x0,y1,z1)),
    tri((-1,0,0), (x0,y0,z0), (x0,y1,z1), (x0,y1,z0)),
]
# +x face, outward normal (+1, 0, 0): go +y then +z
faces += [
    tri((+1,0,0), (x1,y0,z0), (x1,y1,z0), (x1,y1,z1)),
    tri((+1,0,0), (x1,y0,z0), (x1,y1,z1), (x1,y0,z1)),
]
# -y face, outward normal (0, -1, 0): go +x then +z
faces += [
    tri((0,-1,0), (x0,y0,z0), (x1,y0,z0), (x1,y0,z1)),
    tri((0,-1,0), (x0,y0,z0), (x1,y0,z1), (x0,y0,z1)),
]
# +y face, outward normal (0, +1, 0): go +z then +x
faces += [
    tri((0,+1,0), (x0,y1,z0), (x1,y1,z1), (x1,y1,z0)),
    tri((0,+1,0), (x0,y1,z0), (x0,y1,z1), (x1,y1,z1)),
]
# -z face, outward normal (0, 0, -1): go +y then +x
faces += [
    tri((0,0,-1), (x0,y0,z0), (x0,y1,z0), (x1,y1,z0)),
    tri((0,0,-1), (x0,y0,z0), (x1,y1,z0), (x1,y0,z0)),
]
# +z face, outward normal (0, 0, +1): go +x then +y
faces += [
    tri((0,0,+1), (x0,y0,z1), (x1,y1,z1), (x0,y1,z1)),
    tri((0,0,+1), (x0,y0,z1), (x1,y0,z1), (x1,y1,z1)),
]

os.makedirs('constant/triSurface', exist_ok=True)
out = 'constant/triSurface/plate.stl'
with open(out, 'w') as f:
    f.write('solid plate\n')
    for face in faces:
        f.write(face)
    f.write('endsolid plate\n')

print(f'Written {out}  (t = {t:.3e} m)')
