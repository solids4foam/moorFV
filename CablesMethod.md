# METHOD

How to use moorFV beams as restraints in a CFD simulation. Here, two beams are used. Note that beam___ means either beamone or beamtwo. All the steps need to be done for both beamone and beamtwo.

## Step 1:

Copy the  moorFV/beam___ directory to a new separate beamDir/ directory outside the CFD_simulation folder.

## Step 2:

Edit beamDir/ to have correct points, rotation etc... Beam length and radius are specified in beamDir/beam___/constant/beamProperties. The rotation and tranlation are specified in beamDir/beam___/Allrun. 

## Step 3:

Run beamDir/beam___/Allrun scripts. This step generates the beams. They can be viewed in ParaView as a case.foam file will be generated as part of the Allrun script. Note that to view them properly on ParaView, "Warp by Vector" needs to be activated with the vector that the beams are "warped" by set as pointW. 

## Step 4:

Copy 0, constant, system directories from beamDir/beam___ into CFD_simulation/0.orig/beam___, CFD_simulation/constant/beam___, and CFD_simulation/system/beam___.

## Step 5:

Edit CFD_simulation/0.orig/beam___/W and CFD_simulation/0.orig/beam___/Theta. The W files specify the the boundary fields of the beam; both the left and right should have:

type fixedValue;
value uniform (0 0 0);

The Theta files need to be edited to ensure that OpenFOAM is looking in the right place for information regarding the roation. In both the left and the right boundary fields, we should have:

file "$FOAM_CASE/constant/beam___/timeVsMoment";

## Step 6:

Edit CFD_simulation/constant/dynaicMeshDict. Add restraints whilst ensuring that the attachment points are correct.

        restraints
        {
                beamone
        {
                sixDoFRigidBodyMotionFvBeamRestraint finiteVolumeBeam;

                attachmentPatch right;

                anchorPatch left;

                beamRegion      "beamone";

                refAttachmentPt (5.75 2.375 0.8);
        }
            
                beamtwo
        {
                sixDoFRigidBodyMotionFvBeamRestraint finiteVolumeBeam;

                attachmentPatch left;

                anchorPatch     right;

                beamRegion      "beamtwo";

                refAttachmentPt (5.75 2.375 0.8);
        }
	}

Also need to ensure that the following are set corrrectly:

motionSolverLibs    ("libfiniteVolumeBeamMooring");

motionSolver        sixDoFRigidBodyMotionFvBeam;

The title for the coefficients needs to be correct too:

sixDoFRigidBodyMotionFvBeamCoeffs


## Step 7:

Run your CFD_simulation OpenFOAM case and view it in Paraview. Note that it is more straightforward to run the case in serial instead of parallel. Extra steps need to be taken in the CFD_simulations/Allrun file if a parallel run is desired, but it is possible.
