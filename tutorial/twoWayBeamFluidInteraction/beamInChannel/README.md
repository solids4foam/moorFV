# Flexible beam in a confined channel flow

A soft cantilever, clamped to the lower wall of a 2D channel, bends under a parabolic inflow
until it reaches a steady deflection. Reference: Zhang, Liu & Khoo (2012), IJNME,
doi:10.1002/nme.4299, Section 6.4, Case 1.

- `beamTunnel/`: coupled case (`pimpleFoam` + actuator-line beam), Cdn = 26.73
- `beamInitialisation/`: generates the beam mesh and initial beam fields
- `dragCalib/`: steady `simpleFoam` run on a rigid plate, used to obtain Cdn

## Run

```sh
./Allrun
```

This generates the beam in `beamInitialisation/`, copies its mesh and fields into
`beamTunnel/`, and runs `beamTunnel/`. `./Allclean` removes them again. The tip displacement is
extracted with `beamTunnel/bash_tip.sh`.

## Citation

This case is described in:

Taran, A., Bali, S., Tuković, Ž., Pakrashi, V., and Cardiff, P. (2026). Two-way actuator-line
coupling of a geometrically exact finite-volume beam with incompressible and free-surface flow.
arXiv:2609.39443. https://doi.org/10.48550/arXiv.2609.39443
