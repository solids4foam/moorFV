# Using moorFV Beams as Restraints in a CFD Simulation

Here, two beams are used. Note that `beam___` means either beamone or beamtwo. All the steps need to be done for both beamone and beamtwo.

## Step 1:

Copy the  `moorFV/beam___` directory to a new separate `beamDir/` directory outside the CFD_simulation folder.

## Step 2:

Edit `beamDir/` to have correct points, rotation, etc. Beam length, radius and number of segments are specified in the `beams` list of `beamDir/beam___/constant/beamProperties`. Seabed contact, Morison drag and Morison added mass are specified in `beamDir/beam___/constant/beamMomentumContributionProperties`. The initial rotation and translation are specified in `beamDir/beam___/Allrun`.

## Step 3:

Run `beamDir/beam___/Allrun` scripts. This step generates the beams. Then run `beamFoam` in `beamDir/beam___` to let each line settle into its catenary. They can be viewed in ParaView as a case.foam file will be generated as part of the Allrun script. Note that to view them properly on ParaView, **"Warp by Vector"** needs to be activated with the vector that the beams are **"warped"** by set as `pointW`.

## Step 4:

Copy the `0/ref*` fields and the latest time directory from `beamDir/beam___` into `CFD_simulation/0.orig/beam___`. The `ref*` fields hold the initial beam position and are only written to `0/` (by `setInitialPositionBeam`). Copy the `constant/` and `system/` directories into `CFD_simulation/constant/beam___` and `CFD_simulation/system/beam___`. `system/beam___/fvSchemes` must have a `d2dt2Schemes` entry (e.g. `default Euler;`).

In the tutorial cases, `prepareBeams` does Steps 3 to 5 for every beam in `cables_initialization`. The `system/beam___` directories are already in `mooredFloatingObject`.

## Step 5:

Edit `CFD_simulation/0.orig/beam___/W` and `CFD_simulation/0.orig/beam___/Theta`. The `W` indicate the displacement of the beam; both the left and right should be like:
```cpp
type fixedValue;

value uniform (0 0 0);
```

The `Theta` files need to be edited to ensure that OpenFOAM is looking in the right place for information regarding the roation. In both the left and the right boundary fields, it should be like:

```cpp
file "$FOAM_CASE/constant/beam___/timeVsMoment";
```

## Step 6:

Edit `CFD_simulation/constant/dynamicMeshDict`. Add restraints whilst ensuring that the attachment points are correct.
```cpp
restraints
{
    beamone
    {
        sixDoFRigidBodyBeamMotionRestraint finiteVolumeBeam;

        attachmentPatch right;

        anchorPatch left;

        beamRegion      "beamone";

        refAttachmentPt (5.75 2.375 0.8);
    }
    beamtwo
    {
        sixDoFRigidBodyBeamMotionRestraint finiteVolumeBeam;

        attachmentPatch left;

        anchorPatch     right;

        beamRegion      "beamtwo";

        refAttachmentPt (5.75 2.375 0.8);
    }
}
```
Also need to ensure that the following are set corrrectly within :

```cpp
motionSolverLibs    ("libsixDoFRigidBodyBeamMotion");

motionSolver        sixDoFRigidBodyBeamMotion;
```

Please note that the title for the coefficients sub-dict needs to be set as `sixDoFRigidBodyBeamMotionCoeffs`.


## Step 7:

Run your CFD_simulation OpenFOAM case and visualise it in Paraview. Note that it is more straightforward to run the case in serial instead of parallel. Extra steps need to be taken in the `CFD_simulations/Allrun` file if a parallel run is desired, but it is possible.
