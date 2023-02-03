/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "finiteVolumeBeamRestraint.H"
#include "addToRunTimeSelectionTable.H"
#include "sixDoFRigidBodyMotionNew.H"
#include "beamModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sixDoFRigidBodyMotionRestraints
{
    defineTypeNameAndDebug(finiteVolumeBeamRestraint, 0);
    addToRunTimeSelectionTable
    (
        sixDoFRigidBodyMotionRestraintNew,
        finiteVolumeBeamRestraint,
        dictionary
    );
};
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionRestraints::finiteVolumeBeamRestraint::finiteVolumeBeamRestraint
(
    const dictionary& sDoFRBMRDict, const Time& time
)
:
    sixDoFRigidBodyMotionRestraintNew(sDoFRBMRDict, time),
    anchor_(),
    refAttachmentPt_(),
    stiffness_(),
    damping_(),
    restLength_()
{
    read(sDoFRBMRDict);

    Foam::autoPtr<Foam::beamModel> beam =
        Foam::beamModel::New
        (
            const_cast<Time&>(time) , Foam::dynamicFvMesh::defaultRegion
        );

    // beamPtr
    // (
    //     new beamModels::coupledTotalLagNewtonRaphsonBeam
    //     (
    //         const_cast<Time&>(time)
    //     )
    // );
}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionRestraints::finiteVolumeBeamRestraint::~finiteVolumeBeamRestraint()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::sixDoFRigidBodyMotionRestraints::finiteVolumeBeamRestraint::restrain
(
    const sixDoFRigidBodyMotionNew& motion,
    vector& restraintPosition,
    vector& restraintForce,
    vector& restraintMoment
) const
{
    restraintPosition = motion.currentPosition(refAttachmentPt_);

    vector r = restraintPosition - anchor_;

    scalar magR = mag(r);

    // r is now the r unit vector
    r /= (magR + VSMALL);

    vector v = motion.currentVelocity(restraintPosition);

    restraintForce = -stiffness_*(magR - restLength_)*r - damping_*(r & v)*r;

    restraintMoment = vector::zero;

    if (motion.report())
    {
        Info<< " spring length " << magR
            << " force " << restraintForce
            << " moment " << restraintMoment
            << endl;
    }
}


bool Foam::sixDoFRigidBodyMotionRestraints::finiteVolumeBeamRestraint::read
(
    const dictionary& sDoFRBMRDict
)
{
    sixDoFRigidBodyMotionRestraintNew::read(sDoFRBMRDict);

    sDoFRBMRCoeffs_.lookup("anchor") >> anchor_;

    sDoFRBMRCoeffs_.lookup("refAttachmentPt") >> refAttachmentPt_;

    sDoFRBMRCoeffs_.lookup("stiffness") >> stiffness_;

    sDoFRBMRCoeffs_.lookup("damping") >> damping_;

    sDoFRBMRCoeffs_.lookup("restLength") >> restLength_;

    return true;
}


void Foam::sixDoFRigidBodyMotionRestraints::finiteVolumeBeamRestraint::write
(
    Ostream& os
) const
{
    os.writeKeyword("anchor")
        << anchor_ << token::END_STATEMENT << nl;

    os.writeKeyword("refAttachmentPt")
        << refAttachmentPt_ << token::END_STATEMENT << nl;

    os.writeKeyword("stiffness")
        << stiffness_ << token::END_STATEMENT << nl;

    os.writeKeyword("damping")
        << damping_ << token::END_STATEMENT << nl;

    os.writeKeyword("restLength")
        << restLength_ << token::END_STATEMENT << nl;
}

// ************************************************************************* //
