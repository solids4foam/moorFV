---
name: moorFV
description: >
  Use this skill for any task involving the moorFV repository: adding new
  classes, editing motion solvers or restraints, modifying build files,
  writing tutorials, following OpenFOAM runtime selection table patterns,
  or maintaining coding style. Trigger whenever the user mentions moorFV,
  mooring, beamActuatorLine, sixDoFRigidBodyBeamMotion, beamFoam, ALM,
  actuator line, wmake, or requests C++ changes in this codebase.
---

# Codex Guidelines for moorFV

This document defines how automated coding changes should be made in this repository.

## 1) Core Principles

- Make the smallest correct change.
- Preserve existing architecture and naming patterns.
- Prefer consistency with nearby code over introducing new style variants.
- Avoid broad refactors unless explicitly requested.
- Keep compatibility with OpenFOAM.com (v2306 is the primary tested version) in mind.
- State assumptions before implementation when behavior, target OpenFOAM
  variant, or expected numerical result is unclear.
- If multiple interpretations are plausible, ask or present the tradeoff instead
  of silently choosing.
- Push back on requests that would cause unnecessary redesign, broad
  compatibility risk, or speculative behavior.

## 2) Repository Purpose and Scope

moorFV is an OpenFOAM coupling library with two distinct use cases:

1. **Moored floating bodies**: couples the `sixDoFRigidBodyBeamMotion` mesh
   motion solver (6-DoF rigid body dynamics) with beamFoam finite-volume
   Simo-Reissner beam elements acting as mooring lines.

2. **Beam-fluid interaction**: distributes beam-computed actuator-line forces
   to the fluid momentum equation via `beamActuatorLine` (an fvOption using a
   Gaussian projection kernel — not porosity).

The `beamFoam` submodule (located at `src/beamFoam/`, branch `alm-on-main`)
provides the beam solver. moorFV wraps it; do not modify beamFoam files
directly.

## 3) Coding Style Rules

### C++ style

- Follow existing OpenFOAM/moorFV style in surrounding files.
- Use the same indentation, brace style, comment style, and naming conventions
  as local code.
- Keep lines and expressions readable; avoid clever/condensed code.
- Prefer explicit, local, maintainable changes over abstraction-heavy rewrites.
- Add comments only when behavior is non-obvious; do not add redundant comments.

### File/header conventions

- Preserve existing license/header block format in C++ files.
- Keep include ordering consistent with nearby files.
- Do not change copyright headers unless explicitly asked.
- Header guards must match the file name exactly (e.g., `#ifndef beamActuatorLine_H`).
- For files in different namespaces that share a base name, scope the guard
  (e.g., `#ifndef functionObjects_sixDoFRigidBodyBeamMotionState_H`).

### Scripts and docs

- Match existing shell script style in `Allwmake`, `Allwclean`, `Allrun`, etc.
- Keep Markdown concise, practical, and repository-specific.

## 4) OpenFOAM Conventions to Follow

- Respect OpenFOAM runtime type conventions:
  - `TypeName("...")` in class headers.
  - Registration via `addToRunTimeSelectionTable(...)` in source files.
- Preserve dictionary-driven behavior and runtime configurability.
- This repository targets **OpenFOAM.com** style only (not foam-extend or
  OpenFOAM.org); do not add multi-flavor `#ifdef` guards.
- Avoid introducing dependencies that break existing wmake workflows.

## 5) Runtime Selection Tables (How They Work)

moorFV uses OpenFOAM runtime selection tables to instantiate models from
dictionaries at runtime.

Key runtime-selectable types in this repository:

| Type string (in dict) | C++ class | Table |
|---|---|---|
| `sixDoFRigidBodyBeamMotion` | `sixDoFRigidBodyBeamMotion` | `motionSolver` |
| `FvBeamNewmark` | `FvBeamNewmark` | `sixDoFFvBeamSolver` |
| `finiteVolumeBeam` | `finiteVolumeBeam` | `sixDoFRigidBodyBeamMotionRestraint` |
| `beamActuatorLine` | `beamActuatorLine` | `fvOption` |
| `sixDoFRigidBodyBeamMotionState` | `sixDoFRigidBodyBeamMotionState` | `functionObject` |
| Constraint types: `axisFvBeam`, `planeFvBeam`, `pointFvBeam`, `lineFvBeam`, `orientationFvBeam` | respective constraint classes | `sixDoFRigidBodyBeamMotionConstraint` |

Rules when adding new runtime-selectable classes:

- Add `TypeName("...")` in header.
- Add `addToRunTimeSelectionTable(...)` in source.
- Add source file to `src/sixDoFRigidBodyBeamMotion/Make/files`.
- Ensure the dictionary `type` string exactly matches the `TypeName` value.

## 6) Build System

Two libraries are compiled by `src/Allwmake`:

| Directory | Output library |
|---|---|
| `src/sixDoFRigidBodyBeamMotion/` | `$FOAM_USER_LIBBIN/libsixDoFRigidBodyBeamMotion.so` |
| `src/sixDoFRigidBodyBeamMotionState/` | `$FOAM_USER_LIBBIN/libsixDoFRigidBodyBeamMotionState.so` |

Tutorial cases load these in `system/controlDict` under `libs`:

```c++
libs
(
    "libsixDoFRigidBodyBeamMotion.so"
    "libsixDoFRigidBodyBeamMotionState.so"   // if using the state function object
);
```

The motion solver is declared in `constant/dynamicMeshDict`:

```c++
motionSolver        sixDoFRigidBodyBeamMotion;
motionSolverLibs    ("libsixDoFRigidBodyBeamMotion");
```

When adding a source file:
- Add `.C` path to `src/sixDoFRigidBodyBeamMotion/Make/files` (or the state
  library's `Make/files` if it belongs there).
- Never add beamFoam source files to moorFV build lists.

## 7) beamFoam Dependency

- beamFoam is a **git submodule** at `src/beamFoam/` on branch `alm-on-main`.
- moorFV uses beamFoam headers from `$BEAMFOAM_DIR/src/wireBunchingModels/lnInclude`.
- beamFoam must be compiled before moorFV (`cd src/beamFoam && ./Allwmake -j`).
- Do not modify any file under `src/beamFoam/`; it is an external dependency.
- Key beamFoam types used by moorFV:
  - `coupledTotalLagNewtonRaphsonBeam` — the beam solver that computes `W`,
    `almForce`, `fluidCellIDs`, `refW`.
  - `beamModel` — base class for beam models.

## 8) fvOption: beamActuatorLine

`beamActuatorLine` (in `src/sixDoFRigidBodyBeamMotion/fvOptions/beamActuatorLine/`)
distributes per-unit-length beam forces to fluid cells using a 2D Gaussian
kernel. It is **not** a porosity model.

Key implementation points:
- `eta(r)` — 2D Gaussian kernel: `(1/(ε²π)) exp(-(r/ε)²)`
- `calculateS()` — for each fluid cell, finds the closest beam segment (index,
  normalized coordinate `s ∈ [0,1]`, radial distance `r`).
- `applyBeamForce()` — reads `almForce` from the beam mesh, distributes via
  kernel-weighted interpolation, adds `(Sj/fluidRho_)*V[c]` to `eqn.source()`.
- Writes `eta` and `beamActuatorForce` fields at write times.

Required fields on the beam mesh: `refW`, `W`, `almForce`.

**Known limitation**: both `addSup` overloads (incompressible and
density-weighted) call the same `applyBeamForce` which divides by the
hard-coded scalar `fluidRho_`. The density-weighted path (interFoam) should
instead add `Sj*V[c]` without dividing by rho — this is a pending fix.

Dictionary usage:

```c++
beamActuatorLine
{
    type            beamActuatorLine;
    active          yes;

    beamActuatorLineCoeffs
    {
        epsilon     0.05;    // Gaussian kernel half-width [m]
        beamName    beam;    // beam region name
        rho         1000;    // fluid density [kg/m³]
    }
}
```

## 9) Minimise Changes

- Only modify files necessary for the requested task.
- Do not reformat unrelated code.
- Do not rename symbols/files unless required.
- Do not alter behavior outside requested scope.
- Prefer targeted edits over cleanup passes.
- Every changed line should trace directly to the user's request, a required
  build fix, or a test/documentation update needed to verify the change.
- If unrelated dead code or cleanup is noticed, mention it separately; do not
  remove it unless asked.

Before finalizing, verify:

- Build impact is localized.
- No unrelated files changed.

## 10) Change Delivery Format

- Keep changes patch-oriented and reviewable.
- When editing directly in the workspace, summarize changed files and
  verification performed.
- If multiple concerns are required, split into logical commits.

## 11) Practical Workflow

1. Define the concrete success criteria for the request.
2. Read nearby code and follow local patterns.
3. Implement the minimal patch.
4. Update runtime registration and `Make/files` if adding a new class.
5. Compile with `wmake libso sixDoFRigidBodyBeamMotion` (or the full
   `src/Allwmake`) to verify.
6. Report what was verified and what was not.

## 12) What to Avoid

- Large-scale refactors without explicit request.
- API redesigns when a local fix is sufficient.
- New dictionary options, runtime switches, or configurability unless requested.
- Introducing new style conventions inconsistent with repository norms.
- Silent behavioral changes not documented in the patch summary.
- Modifying files inside `src/beamFoam/`.

## 13) Reference Files (Commonly Relevant)

- Motion solver: `src/sixDoFRigidBodyBeamMotion/sixDoFRigidBodyBeamMotionSolver/sixDoFRigidBodyBeamMotionSolver.H`
- Main 6-DoF class: `src/sixDoFRigidBodyBeamMotion/sixDoFRigidBodyBeamMotion/sixDoFRigidBodyBeamMotion.H`
- Beam restraint: `src/sixDoFRigidBodyBeamMotion/sixDoFRigidBodyBeamMotion/restraints/finiteVolumeBeamRestraint/finiteVolumeBeam.H`
- Newmark ODE solver: `src/sixDoFRigidBodyBeamMotion/sixDoFFvBeamSolvers/FvBeamNewmark/FvBeamNewmark.H`
- Actuator line fvOption: `src/sixDoFRigidBodyBeamMotion/fvOptions/beamActuatorLine/beamActuatorLine.H`
- State function object: `src/sixDoFRigidBodyBeamMotionState/sixDoFRigidBodyBeamMotionState.H`
- Point patch BC: `src/sixDoFRigidBodyBeamMotion/FvBeamPointPatchFields/derived/finiteVolumeSixDoFRigidBodyDisplacement/finiteVolumeSixDoFRigidBodyDisplacementPointPatchVectorField.H`
- Build lists: `src/sixDoFRigidBodyBeamMotion/Make/files`, `src/sixDoFRigidBodyBeamMotionState/Make/files`
- Build scripts: `Allwmake`, `src/Allwmake`, `src/Allwclean`
- Tutorial mooring case: `tutorial/mooringCase-H12T20/mooredFloatingObject/`
- Tutorial beam-fluid case: `tutorial/beamFluidInteractionCases/pitzDailyWithBeam/`
