# Leaflet in a sinusoidal channel flow

A thin flexible leaflet, clamped at the channel mid-plane, oscillates in a sinusoidal channel
flow. Reference: Zhang & Gay (2007), J. Fluids Struct. 23, 839-857, Section 3.2.

| Case | Re | St | Period T (s) |
|---|---|---|---|
| `C1` | 10 | 1.0 | 1 |
| `C2` | 1 | 1.0 | 1 |
| `C3` | 1 | 0.5 | 2 |
| `C4` | 10 | 0.5 | 2 |

`dragCalib/` holds the steady `simpleFoam` runs used to obtain Cdn: `topView` (cross-section in
free stream) and `sideView-Re10`, `sideView-Re1` (rigid plate in the channel).

## Run

```sh
cd C1
./Allrun
```

The tip displacement is extracted with `bash_tip.sh`.
