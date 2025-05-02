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
    patchID_(-1),
    anchorPatchID_(-1),
    initialW_(vector::zero),
    storeInitialW_(true),
    initialQ_(vector::zero),
    forceFilePtr_(),
    anchorForceFilePtr_()
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

    }

    // Write to the state dictionary to allow restarts
    state_.add("refAttachmentPt", refAttachmentPt_);
    state_.add("attachmentPatch", attachmentPatch_);
    state_.add("anchorPatch", anchorPatch_);
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
    const sixDoFRigidBodyMotionFvBeam& motion,
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
        if (beam.solutionW().boundaryField()[patchID_].size() != 0)
        {
            initialW_ = beam.solutionW().boundaryField()[patchID_][0];
        }

        // Update the state
        state_.add("initialW", initialW_);
        state_.add("storeInitialW", storeInitialW_);
    }

    restraintPosition = motion.transform(refAttachmentPt_);

    const vector attachmentDisp = restraintPosition - refAttachmentPt_;

    volVectorField& W = beam.solutionW();

    // consider relaxing the displacement...
    W.boundaryFieldRef()[patchID_] == attachmentDisp + initialW_;
    beam.evolve();

    beam.updateTotalFields();

    vector localAnchorForce(vector::zero);
    vector localAttachmentForce(vector::zero);
    // Lookup the surface forces
    const surfaceVectorField& Q =
        beam.mesh().lookupObject<surfaceVectorField>("Q");

    if (Q.boundaryField()[patchID_].size() == 1)
    {
        localAnchorForce = Q.boundaryField()[anchorPatchID_][0];
        localAttachmentForce = Q.boundaryField()[patchID_][0];
        // FatalError
        //     << "Q.boundaryField()[patchID_].size() != 1: "
        //     << "this needs to be parallelised" << abort(FatalError);
    }
    reduce(localAttachmentForce, sumOp<vector>());
    reduce(localAnchorForce,     sumOp<vector>());

    vector anchorForce     = localAnchorForce;
    vector attachmentForce = localAttachmentForce;

    restraintForce = -attachmentForce;
    restraintMoment = vector::zero;


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
