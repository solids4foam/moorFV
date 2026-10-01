/*---------------------------------------------------------------------------*\

  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/
#include "finiteVolumeBeam.H"
#include "addToRunTimeSelectionTable.H"
#include "sixDoFRigidBodyMotionFvBeam.H"
#include "Time.H"
#include "fvMesh.H"
#include "OFstream.H"
#include "quaternion.H"
#include "PstreamReduceOps.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sixDoFRigidBodyMotionFvBeamRestraints
{
    defineTypeNameAndDebug(finiteVolumeBeam, 0);

    addToRunTimeSelectionTable
    (
        sixDoFRigidBodyMotionFvBeamRestraint,
        finiteVolumeBeam,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionFvBeamRestraints::finiteVolumeBeam::finiteVolumeBeam
(
    const word& name,
    const dictionary& sDoFRBMRDict,
    const Time& time
)
:
    sixDoFRigidBodyMotionFvBeamRestraint(name, sDoFRBMRDict,time),
    beam_
    (
        beamModel::New
        (
            const_cast<Time&>(time), word(sDoFRBMRCoeffs_.lookup("beamRegion"))
        )
    ),
    state_
    (
        IOobject
        (
            "stateFiniteVolumeBeam" + name,
            time.timeName(),
            "uniform",
            time,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        dictionary()
    ),
    refAttachmentPt_(),
    attachmentPatch_(),
    anchorPatch_(),
    beamAttachmentKinematics_("prescribedByRigidBody"),
    patchID_(-1),
    anchorPatchID_(-1),
    initialW_(vector::zero),
    storeInitialW_(true),
    initialQ_(vector::zero),
    forceFilePtr_(),
    anchorForceFilePtr_(),
    attachmentForceFilePtr_(),
    displacementFilePtr_(),
    monolithicInitialised_(false)
{
    if (debug)
    {
        InfoIn("finiteVolumeBeam::finiteVolumeBeam(...)")
            << "Creating finiteVolumeBeam " << name << endl;
    }

    // Read settings
    read(sDoFRBMRDict);

    // If needed, update patch indices
    if (patchID_ == -1)
    {
        patchID_ =
            beam_->mesh().boundaryMesh().findPatchID
            (
                attachmentPatch_
            );
    }

    if (anchorPatchID_ == -1)
    {
        anchorPatchID_ =
            beam_->mesh().boundaryMesh().findPatchID
            (
                anchorPatch_
            );
    }

    // Create force file
    if (Pstream::master())
    {
        fileName historyDir;

        word startTimeName =
            time.timeName(time.startTime().value());

        if (Pstream::parRun())
        {
            // Put in undecomposed case (Note: gives problems for
            // distributed data running)
            historyDir = time.path()/".."/"postProcessing"/startTimeName;
        }
        else
        {
            historyDir = time.path()/"postProcessing"/startTimeName;
        }

        // Create directory if does not exist.
        mkDir(historyDir);

        // Open new file at start up
        forceFilePtr_.reset
        (
            new OFstream
            (
                historyDir/"force" + name + ".dat"
            )
        );
        anchorForceFilePtr_.reset
        (
            new OFstream
            (
                historyDir/"anchorForce" + name + ".dat"
            )
        );
        attachmentForceFilePtr_.reset
        (
            new OFstream
            (
                historyDir/"attachmentForce" + name + ".dat"
            )
        );
        
        displacementFilePtr_.reset
        (
            new OFstream
            (
                historyDir/"displacement" + name + ".dat"
            )
        );


        // Add headers to output data
        if (forceFilePtr_.valid())
        {
            forceFilePtr_()
                << "# Time"
                << " " << "forceX"
                << " " << "forceY"
                << " " << "forceZ"
                << endl;
        }
        if (anchorForceFilePtr_.valid())
        {
            anchorForceFilePtr_()
                << "# Time"
                << " " << "anchorForceX"
                << " " << "anchorForceY"
                << " " << "anchorForceZ"
                << endl;
        }
        if (attachmentForceFilePtr_.valid())
        {
            attachmentForceFilePtr_()
                << "# Time"
                << " " << "attachmentForceX"
                << " " << "attachmentForceY"
                << " " << "attachmentForceZ"
                << endl;
        }
        if (displacementFilePtr_.valid())
        {
            displacementFilePtr_()
                << "# Time"
                << " " << "dispX"
                << " " << "dispY"
                << " " << "dispZ"
                << endl;
        }

    }

    // Write to the state dictionary to allow restarts
    state_.add("refAttachmentPt", refAttachmentPt_);
    state_.add("attachmentPatch", attachmentPatch_);
    state_.add("anchorPatch", anchorPatch_);
    state_.add("beamAttachmentKinematics", beamAttachmentKinematics_);
    state_.add("patchID", patchID_);
    state_.add("anchorPatchID", anchorPatchID_);
    state_.add("initialQ", initialQ_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionFvBeamRestraints::finiteVolumeBeam::
~finiteVolumeBeam()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


void Foam::sixDoFRigidBodyMotionFvBeamRestraints::finiteVolumeBeam::restrain
(
    sixDoFRigidBodyMotionFvBeam& motion,
    vector& restraintPosition,
    vector& restraintForce,
    vector& restraintMoment
) const
{
    if (debug)
    {
        InfoIn("finiteVolumeBeam::restrain(...)")
            << endl;
    }

    // Take a reference to the beam model
    beamModel& beam = beam_();
    const scalar t = motion.time().timeOutputValue();
    if (storeInitialW_)
    {
        Info<< "Storing initial W" << endl;
        storeInitialW_ = false;
        if (beam.solutionW().boundaryField()[patchID_].size() == 0)
        {
            FatalError
                << "W.boundaryField()[patchID_].size() == 0!"
                << abort(FatalError);
        }
        initialW_ = beam.solutionW().boundaryField()[patchID_][0];

        // Update the state
        state_.add("initialW", initialW_);
        state_.add("storeInitialW", storeInitialW_);
    }

    restraintPosition = motion.transform(refAttachmentPt_);

    const vector attachmentDisp = restraintPosition - refAttachmentPt_;

    volVectorField& W = beam.solutionW();

    if (beamAttachmentKinematics_ == "prescribedByRigidBody")
    {
        // consider relaxing the displacement...
        W.boundaryFieldRef()[patchID_] == attachmentDisp + initialW_;
    }
    else if (beamAttachmentKinematics_ == "solvedByBlockEigen")
    {
        Info<< "finiteVolumeBeam: attachment kinematics owned by BlockEigen"
            << endl;
    }
    else
    {
        FatalErrorInFunction
            << "Unknown beamAttachmentKinematics "
            << beamAttachmentKinematics_ << nl
            << "Valid options are prescribedByRigidBody and solvedByBlockEigen"
            << abort(FatalError);
    }

    //---------------------------------------------------------------------
    // Colm: include all rigid body data- not just mooring attachment points
    //---------------------------------------------------------------------
    RigidBodyStepData rbData;

    const sixDoFRigidBodyMotionFvBeamState& state = motion.state();
    const sixDoFRigidBodyMotionFvBeamState& state0 = motion.state0();
    const sixDoFRigidBodyMotionFvBeam& constMotion = motion;

    rbData.previous.displacement =
        state0.centreOfRotation() - constMotion.initialCentreOfRotation();

    rbData.current.displacement =
        state.centreOfRotation() - constMotion.initialCentreOfRotation();

    rbData.previous.orientation = state0.Q();
    rbData.current.orientation = state.Q();

    rbData.previous.velocity = state0.v();
    rbData.current.velocity = state.v();

    rbData.previous.angularMomentum = state0.pi();
    rbData.current.angularMomentum = state.pi();

    rbData.previous.acceleration = state0.a();
    rbData.current.acceleration = state.a();

    rbData.previous.torque = state0.tau();
    rbData.current.torque = state.tau();

    rbData.mass = motion.mass();
    rbData.momentOfInertia = motion.momentOfInertia();
    rbData.initialCentreOfRotation = constMotion.initialCentreOfRotation();
    rbData.centreOfRotation = motion.centreOfRotation();
    rbData.attachmentPoint = restraintPosition;
    rbData.solveAttachmentKinematicsInBlockEigen =
        beamAttachmentKinematics_ == "solvedByBlockEigen";

    //    beamModels::coupledTotalLagNewtonRaphsonBeam& coupledBeam =
    //    refCast<beamModels::coupledTotalLagNewtonRaphsonBeam>(beam);

    // coupledBeam.setRigidBodyStepData(rbData);
      
    beam.setRigidBodyStepData(rbData);
 
//---------------------------------------------------------------------    

    // Call beamFoam- Colm
    beam.evolve();

    const Switch blockEigenAuthoritativeMotion
    (
        sDoFRBMRCoeffs_.getOrDefault<Switch>
        (
            "blockEigenAuthoritativeMotion",
            false
        )
    );

    RigidBodySolution rigidBodySolution;
    bool blockEigenMotionApplied = false;

    if (blockEigenAuthoritativeMotion)
    {
        if (beam.getRigidBodySolution(rigidBodySolution))
        {
            motion.applyBlockEigenSolution(rigidBodySolution);
            restraintPosition = motion.transform(refAttachmentPt_);
            blockEigenMotionApplied = true;

            Info<< "finiteVolumeBeam using authoritative BlockEigen "
                << "rigid-body solution" << nl
                << "    centreOfRotation = "
                << motion.centreOfRotation() << nl
                << "    velocity = " << motion.v() << nl
                << "    acceleration = " << motion.state().a()
                << endl;
        }
        else
        {
            WarningInFunction
                << "blockEigenAuthoritativeMotion requested, but the beam "
                << "model did not provide a valid rigid-body solution"
                << endl;
        }
    }

    beam.updateTotalFields();

    // Lookup the surface forces
    const surfaceVectorField& Q =
        beam.mesh().lookupObject<surfaceVectorField>("Q");

    if (Q.boundaryField()[patchID_].size() != 1)
    {
        FatalError
            << "Q.boundaryField()[patchID_].size() != 1: "
            << "this needs to be parallelised" << abort(FatalError);
    }

    const vector attachmentForce = Q.boundaryField()[patchID_][0];
    const vector anchorForce = Q.boundaryField()[anchorPatchID_][0];
    const vector anchorDisplacement = W.boundaryField()[anchorPatchID_][0];  // first cell of patch

    const vector beamReactionForce = -attachmentForce;
    restraintForce = beamReactionForce;
    // relax force
    // store and lookup alpha from dict
    //restraintForce = alpha*restraintForce + (1 - alpha)*restraintForcePrevious;
    restraintMoment = vector::zero;

    if (blockEigenMotionApplied)
    {
        restraintForce = vector::zero;
        restraintMoment = vector::zero;
    }

    if (forceFilePtr_.valid())
    {
        forceFilePtr_()
            << t
            << " " << restraintForce.x()
            << " " << restraintForce.y()
            << " " << restraintForce.z()
            << endl;
    }
    if (anchorForceFilePtr_.valid())
    {
        anchorForceFilePtr_()
            << t
            << " " << anchorForce.x()
            << " " << anchorForce.y()
            << " " << anchorForce.z()
            << endl;
    }
    if (attachmentForceFilePtr_.valid())
    {
        attachmentForceFilePtr_()
            << t
            << " " << attachmentForce.x()
            << " " << attachmentForce.y()
            << " " << attachmentForce.z()
            << endl;
    }
    if (displacementFilePtr_.valid())
    {
        displacementFilePtr_()
            << t
            << " " << anchorDisplacement.x()
            << " " << anchorDisplacement.y()
            << " " << anchorDisplacement.z()
            << endl;
    }

}


Foam::RigidBodyEndState
Foam::sixDoFRigidBodyMotionFvBeamRestraints::finiteVolumeBeam::solveMonolithic
(
    const sixDoFRigidBodyMotionFvBeam& motion,
    const vector& force,
    const vector& moment,
    const dictionary& solverDict
) const
{
    beamModel& beam = beam_();
    const scalar t = motion.time().timeOutputValue();

    if (!monolithicInitialised_)
    {
        monolithicInitialised_ = true;

        const point centreOfMass = motion.centreOfMass();

        if (mag(motion.centreOfRotation() - centreOfMass) > SMALL)
        {
            FatalErrorInFunction
                << "beamFoamCoupled needs the centre of rotation at the "
                << "centre of mass; got centreOfRotation "
                << motion.centreOfRotation() << ", centreOfMass "
                << centreOfMass << abort(FatalError);
        }

        if (mag(centreOfMass - motion.initialCentreOfMass()) > SMALL)
        {
            FatalErrorInFunction
                << "beamFoamCoupled must start from the initial body state; "
                << "restarts are not supported yet" << abort(FatalError);
        }

        // The body in beamFoam: moorFV passes the weight in the applied
        // force, so beamFoam applies no gravity of its own
        dictionary bodyDict;
        bodyDict.add("patch", attachmentPatch_);
        bodyDict.add("mass", motion.mass());
        bodyDict.add("momentOfInertia", motion.momentOfInertia());
        bodyDict.add("orientation", motion.orientation());
        bodyDict.add("centreOfMass", centreOfMass);
        bodyDict.add("velocity", motion.v());
        bodyDict.add("angularVelocity", motion.omega());
        bodyDict.add("newmarkBeta", solverDict.getOrDefault<scalar>("beta", 0.25));
        bodyDict.add("newmarkGamma", solverDict.getOrDefault<scalar>("gamma", 0.5));
        bodyDict.add
        (
            "linearDamping",
            solverDict.getOrDefault<scalar>("linearDamping", 0)
        );
        bodyDict.add
        (
            "angularDamping",
            solverDict.getOrDefault<scalar>("angularDamping", 0)
        );
        bodyDict.add
        (
            "jacobianCheck",
            solverDict.getOrDefault<label>("jacobianCheck", 0)
        );

        beam.initialiseRigidBodyEnd(bodyDict, vector::zero);

        // The attachment point beamFoam uses is the beam end; the moorFV
        // attachment point should coincide with it
        const point beamEnd =
            beam.mesh().Cf().boundaryField()[patchID_][0]
          + beam.solutionW().boundaryField()[patchID_][0];

        if (mag(beamEnd - refAttachmentPt_) > 1e-6*mag(beamEnd))
        {
            WarningInFunction
                << "refAttachmentPt " << refAttachmentPt_
                << " differs from the beam end " << beamEnd
                << "; beamFoam attaches the body at the beam end" << endl;
        }
    }

    beam.setRigidBodyEndExternalLoad(force, moment);

    beam.evolve();

    beam.updateTotalFields();

    RigidBodyEndState state;
    if (!beam.getRigidBodyEndState(state))
    {
        FatalErrorInFunction
            << "The beam model returned no rigidBodyEnd state"
            << abort(FatalError);
    }

    // Histories, as for the partitioned restraint: forcebeam is the force
    // from the beam on the body
    const surfaceVectorField& Q =
        beam.mesh().lookupObject<surfaceVectorField>("Q");
    const volVectorField& W = beam.solutionW();

    const vector attachmentForce = Q.boundaryField()[patchID_][0];
    const vector anchorForce = Q.boundaryField()[anchorPatchID_][0];
    const vector anchorDisplacement = W.boundaryField()[anchorPatchID_][0];

    auto writeVector = [t](autoPtr<OFstream>& os, const vector& v)
    {
        if (os.valid())
        {
            os() << t << " " << v.x() << " " << v.y() << " " << v.z() << endl;
        }
    };

    writeVector(forceFilePtr_, state.beamForce);
    writeVector(anchorForceFilePtr_, anchorForce);
    writeVector(attachmentForceFilePtr_, attachmentForce);
    writeVector(displacementFilePtr_, anchorDisplacement);

    return state;
}


bool Foam::sixDoFRigidBodyMotionFvBeamRestraints::finiteVolumeBeam::read
(
    const dictionary& sDoFRBMRDict
)
{
    sixDoFRigidBodyMotionFvBeamRestraint::read(sDoFRBMRDict);

    sDoFRBMRCoeffs_.readEntry("refAttachmentPt", refAttachmentPt_);
    sDoFRBMRCoeffs_.readEntry("attachmentPatch", attachmentPatch_);
    sDoFRBMRCoeffs_.readEntry("anchorPatch", anchorPatch_);
    beamAttachmentKinematics_ =
        sDoFRBMRCoeffs_.getOrDefault<word>
        (
            "beamAttachmentKinematics",
            "prescribedByRigidBody"
        );
    sDoFRBMRCoeffs_.readIfPresent("patchID", patchID_);
    sDoFRBMRCoeffs_.readIfPresent("anchorPatchID", anchorPatchID_);
    sDoFRBMRCoeffs_.readIfPresent("initialW", initialW_);
    sDoFRBMRCoeffs_.readIfPresent("storeInitialW", storeInitialW_);
    sDoFRBMRCoeffs_.readIfPresent("initialQ", initialQ_);

    // Read from state dictionary
    // Note: this typically only contains data when restarting the simulation
    state_.readIfPresent("refAttachmentPt", refAttachmentPt_);
    state_.readIfPresent("attachmentPatch", attachmentPatch_);
    state_.readIfPresent("anchorPatch", anchorPatch_);
    state_.readIfPresent
    (
        "beamAttachmentKinematics",
        beamAttachmentKinematics_
    );
    state_.readIfPresent("patchID", patchID_);
    state_.readIfPresent("anchorPatchID", anchorPatchID_);
    state_.readIfPresent("initialW", initialW_);
    state_.readIfPresent("storeInitialW", storeInitialW_);
    state_.readIfPresent("initialQ", initialQ_);

    return true;
}


void Foam::sixDoFRigidBodyMotionFvBeamRestraints::finiteVolumeBeam::write
(
    Ostream& os
) const
{
    os.writeEntry("refAttachmentPt", refAttachmentPt_);
    os.writeEntry("attachmentPatch", attachmentPatch_);
    os.writeEntry("anchorPatch", anchorPatch_);
    os.writeEntry("patchID", patchID_);
    os.writeEntry("anchorPatchID", anchorPatchID_);
    os.writeEntry("initialW", initialW_);
    os.writeEntry("storeInitialW", storeInitialW_);
    os.writeEntry("initialQ", initialQ_);

    // Write to the state dictionary to allow restarts
    state_.add("refAttachmentPt", refAttachmentPt_);
    state_.add("attachmentPatch", attachmentPatch_);
    state_.add("anchorPatch", anchorPatch_);
    state_.add("patchID", patchID_);
    state_.add("anchorPatchID", anchorPatchID_);
    state_.add("initialW", initialW_);
    state_.add("storeInitialW", storeInitialW_);
    state_.add("initialQ", initialQ_);

    // Info<< nl << "state is " << state_ << endl;
}


// ************************************************************************* //
