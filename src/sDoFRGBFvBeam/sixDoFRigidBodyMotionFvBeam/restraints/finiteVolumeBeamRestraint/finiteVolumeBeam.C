/*---------------------------------------------------------------------------*\
  =========                 |
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
    refAttachmentPt_(),
    attachmentPatch_(),
    patchID_(-1),
    initialW_(vector::zero),
    storeInitialW_(true),
    initialQ_(vector::zero),
    storeInitialQ_(true),
    forceFilePtr_()
{
    if (debug)
    {
        InfoIn("finiteVolumeBeam::finiteVolumeBeam(...)")
            << "Creating finiteVolumeBeam " << name << endl;
    }

    read(sDoFRBMRDict);

    patchID_ =
        beam_->mesh().boundaryMesh().findPatchID
        (
            attachmentPatch_
	    );

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
    }
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
        if (beam.solutionW().boundaryField()[patchID_].size() == 0)
        {
            FatalError
                << "W.boundaryField()[patchID_].size() == 0!"
                << abort(FatalError);
        }
        initialW_ = beam.solutionW().boundaryField()[patchID_][0];
    }

    restraintPosition = motion.transform(refAttachmentPt_);

    const vector attachmentDisp = restraintPosition - refAttachmentPt_;

    volVectorField& W = beam.solutionW();

    // consider relaxing the displacement...
    W.boundaryFieldRef()[patchID_] == attachmentDisp + initialW_;

    beam.evolve();

    beam.updateTotalFields();
    //beam.writeFields();
    Info<< "patchID_ :" << patchID_ << endl;

    Info<< "attachment patch :" << attachmentPatch_ << endl;

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

    restraintForce = -attachmentForce;
    // relax force
    // store and lookup alpha from dict
    //restraintForce = alpha*restraintForce + (1 - alpha)*restraintForcePrevious;
    restraintMoment = vector::zero;

    Info<< "attachment force = " << attachmentForce << endl;

    if (motion.report())
    {
        Info<< t << " " << restraintForce.x()   //2
            << " " << restraintForce.y()        //3
            << " " << restraintForce.z()        //4
            << endl;
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

}


bool Foam::sixDoFRigidBodyMotionFvBeamRestraints::finiteVolumeBeam::read
(
    const dictionary& sDoFRBMRDict
)
{
    sixDoFRigidBodyMotionFvBeamRestraint::read(sDoFRBMRDict);

    sDoFRBMRCoeffs_.readEntry("refAttachmentPt", refAttachmentPt_);

    sDoFRBMRCoeffs_.readEntry("attachmentPatch", attachmentPatch_);

    return true;
}


void Foam::sixDoFRigidBodyMotionFvBeamRestraints::finiteVolumeBeam::write
(
    Ostream& os
) const
{}


// ************************************************************************* //
