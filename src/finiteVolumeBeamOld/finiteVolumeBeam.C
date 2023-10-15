/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    anchor_(),
    refAttachmentPt_(),
    stiffness_(),
    damping_(),
    restLength_()
    //beam_()
{
    read(sDoFRBMRDict);


    // Make the beam solver object
    // Foam::autoPtr<Foam::beamModel> beam =
    //    Foam::beamModel::New(runTime, Foam::dynamicFvMesh::defaultRegion);
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

    // Calculate new attachment point position
    restraintPosition = motion.transform(refAttachmentPt_);

    // Update the beam solver boundary conditions: set the displacement on the
    // attachment point boundary condition
    // Lookup the W or DW (displacement) field and then find the attachment end
    // boundary condition. Then update the displacement on this boundary
    // condition. 
    
    // Solve the beam equation
    // beam.evolve();

    // Extract the beam end forces and moments







    
    restraintPosition = motion.transform(refAttachmentPt_);

    vector r = restraintPosition - anchor_;

    scalar magR = mag(r);
    r /= (magR + VSMALL);

    vector v = motion.velocity(restraintPosition);

    restraintForce = -stiffness_*(magR - restLength_)*r - damping_*(r & v)*r;

    restraintMoment = Zero;

    if (motion.report())
    {
        Info<< " attachmentPt - anchor " << r*magR
            << " spring length " << magR
            << " force " << restraintForce
            << endl;
    }
}


bool Foam::sixDoFRigidBodyMotionRestraints::finiteVolumeBeam::read
(
    const dictionary& sDoFRBMRDict
)
{
    sixDoFRigidBodyMotionRestraint::read(sDoFRBMRDict);

    sDoFRBMRCoeffs_.readEntry("anchor", anchor_);
    sDoFRBMRCoeffs_.readEntry("refAttachmentPt", refAttachmentPt_);
    sDoFRBMRCoeffs_.readEntry("stiffness", stiffness_);
    sDoFRBMRCoeffs_.readEntry("damping", damping_);
    sDoFRBMRCoeffs_.readEntry("restLength", restLength_);

    return true;
}


void Foam::sixDoFRigidBodyMotionRestraints::finiteVolumeBeam::write
(
    Ostream& os
) const
{
    os.writeEntry("anchor", anchor_);
    os.writeEntry("refAttachmentPt", refAttachmentPt_);
    os.writeEntry("stiffness", stiffness_);
    os.writeEntry("damping", damping_);
    os.writeEntry("restLength", restLength_);
}

// ************************************************************************* //
