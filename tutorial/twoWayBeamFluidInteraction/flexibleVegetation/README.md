# Flexible vegetation stem in regular waves

A flexible blade (Mimic 3, length 0.30 m) mounted on the plateau of a wave flume sways under
regular waves (H = 0.04 m, T = 2.0 s), solved with `interFoam`. Reference: Jacobsen et al.
(2019), J. Fluid Mech. 880, 1036-1069.

Other wave conditions are run by changing `constant/waveProperties` and `Cdn` in
`constant/beam_0/beamMomentumContributionProperties`.

## Run

```sh
./Allrun
```

This generates the beam in `beamInitialisation/`, copies its mesh and fields into
`mimic3-H040T20-cd4/`, and runs it. `./Allclean` removes them again.

The case runs on 32 ranks. If you change this, update both `system/decomposeParDict` and
`system/beam_0/decomposeParDict`. The whole beam is kept on processor 0.

## Citation

This case is described in:

Taran, A., Bali, S., Tuković, Ž., Pakrashi, V., and Cardiff, P. (2026). Two-way actuator-line
coupling of a geometrically exact finite-volume beam with incompressible and free-surface flow.
arXiv:2609.39443. https://doi.org/10.48550/arXiv.2609.39443
