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
#include "sixDoFRigidBodyMotion.H"
#include "Time.H"
#include "fvMesh.H"
#include "OFstream.H"
#include "quaternion.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sixDoFRigidBodyMotionRestraints
{
    defineTypeNameAndDebug(finiteVolumeBeam, 0);

    addToRunTimeSelectionTable
    (
        sixDoFRigidBodyMotionRestraint,
        finiteVolumeBeam,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionRestraints::finiteVolumeBeam::finiteVolumeBeam
(
    const word& name,
    const dictionary& sDoFRBMRDict
)
:
	sixDoFRigidBodyMotionRestraint(name, sDoFRBMRDict),
	beamPtr_(),
	refAttachmentPt_(),
	attachmentPatch_(),
	beamRegion_(),
	patchID_(-1)

{
    read(sDoFRBMRDict);
    patchID_ =
        beamPtr_->mesh().boundaryMesh().findPatchID
        (
            attachmentPatch_
	);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionRestraints::finiteVolumeBeam::~finiteVolumeBeam()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


void Foam::sixDoFRigidBodyMotionRestraints::finiteVolumeBeam::restrain
(
    const sixDoFRigidBodyMotion& motion,
    vector& restraintPosition,
    vector& restraintForce,
    vector& restraintMoment
) const
{
   // We need to update this function to call the FV beam solver

    // Create the beam model if it does not exist
    if (!beamPtr_.valid())
    {
        beamPtr_ =
            Foam::beamModel::New
            (
                const_cast<Time&>(motion.time()),
                Foam::dynamicFvMesh::defaultRegion
		// construct with word(beamRegion_)
           );
    }

    // Take a reference to the beam model
    beamModel& beam = beamPtr_();

    restraintPosition = motion.transform(refAttachmentPt_);

    const vector attachmentDisp = restraintPosition - refAttachmentPt_;

    volVectorField& W = beam.solutionW();

    //W.boundaryField()[patchID_] = attachmentDisp;

    beam.evolve();

    beam.updateTotalFields();

    beam.writeFields();

    const surfaceVectorField& Q =
	    beam.mesh().lookupObject<surfaceVectorField>("Q");

    if (Q.boundaryField()[patchID_].size() != 1)
    {
        FatalError
            << "Q.boundaryField()[patchID_].size() != 1: "
            << "this needs to be parallelised" << abort(FatalError);
    }

    const vector attachmentForce = Q.boundaryField()[patchID_][0];

    restraintForce = attachmentForce;

    restraintMoment = vector::zero;

    Info<< "attachment force = " << attachmentForce << endl;

    if (motion.report())
    {
        Info<< " force " << restraintForce
            << " moment " << restraintMoment
            << endl;
    }


}


bool Foam::sixDoFRigidBodyMotionRestraints::finiteVolumeBeam::read
(
    const dictionary& sDoFRBMRDict
)
{
    sixDoFRigidBodyMotionRestraint::read(sDoFRBMRDict);
    sDoFRBMRCoeffs_.readEntry("refAttachmentPt", refAttachmentPt_);

    sDoFRBMRCoeffs_.readEntry("attachmentPatch", attachmentPatch_);

    sDoFRBMRCoeffs_.readEntry("beamRegion", beamRegion_);

    return true;
}


void Foam::sixDoFRigidBodyMotionRestraints::finiteVolumeBeam::write
(
    Ostream& os
) const
{
}


// ************************************************************************* //
