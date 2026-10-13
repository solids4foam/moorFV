# Flexible vegetation stem in regular waves

A flexible blade (Mimic 3, length 0.30 m) mounted on the plateau of a wave flume sways under
regular waves (H = 0.04 m, T = 2.0 s), solved with `interFoam`. Reference: Jacobsen et al.
(2019), J. Fluid Mech. 880, 1036-1069.

Other wave conditions are run by changing `constant/waveProperties` and `Cdn` in
`constant/beam_0/beamMomentumContributionProperties`.

## Run

```sh
cd mimic3-H040T20-cd4
./Allrun.decompose
```

The case runs on 32 ranks. If you change this, update both `system/decomposeParDict` and
`system/beam_0/decomposeParDict`, and regenerate `constant/cellDecomposition.beams`.
