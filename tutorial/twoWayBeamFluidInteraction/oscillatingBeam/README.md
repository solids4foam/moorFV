# Leaflet in a sinusoidal channel flow

A thin flexible leaflet, clamped at the channel mid-plane, oscillates in a sinusoidal channel
flow. Reference: Zhang & Gay (2007), J. Fluids Struct. 23, 839-857, Section 3.2.

| Case | Re | St | Period T (s) |
|---|---|---|---|
| `C1` | 10 | 1.0 | 1 |
| `C2` | 1 | 1.0 | 1 |
| `C3` | 1 | 0.5 | 2 |
| `C4` | 10 | 0.5 | 2 |

`beamInitialisation/` generates the beam mesh and initial beam fields, shared by all four
cases. `dragCalib/` holds the steady `simpleFoam` runs used to obtain Cdn: `topView`
(cross-section in free stream) and `sideView-Re10`, `sideView-Re1` (rigid plate in the channel).

## Run

```sh
./Allrun C1
```

This generates the beam, copies its mesh and fields into `C1/`, and runs `C1/`. Without
arguments, all four cases are run. `./Allclean` cleans all cases. The tip displacement is
extracted with `bash_tip.sh`.

## Citation

This case is described in:

Taran, A., Bali, S., Tuković, Ž., Pakrashi, V., and Cardiff, P. (2026). Two-way actuator-line
coupling of a geometrically exact finite-volume beam with incompressible and free-surface flow.
arXiv:2609.39443. https://doi.org/10.48550/arXiv.2609.39443
