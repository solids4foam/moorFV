# Flexible beam in a confined channel flow

A soft cantilever, clamped to the lower wall of a 2D channel, bends under a parabolic inflow
until it reaches a steady deflection. Reference: Zhang, Liu & Khoo (2012), IJNME,
doi:10.1002/nme.4299, Section 6.4, Case 1.

- `beamTunnel/`: coupled case (`pimpleFoam` + actuator-line beam), Cdn = 26.73
- `beamInitialisation/`: generates the beam mesh
- `dragCalib/`: steady `simpleFoam` run on a rigid plate, used to obtain Cdn

## Run

```sh
cd beamTunnel
./Allrun
```

To regenerate the beam mesh first, run `./Allclean` and `./Allrun` in this directory, then
`./Allrun` in `beamTunnel/`. The tip displacement is extracted with `beamTunnel/bash_tip.sh`.
