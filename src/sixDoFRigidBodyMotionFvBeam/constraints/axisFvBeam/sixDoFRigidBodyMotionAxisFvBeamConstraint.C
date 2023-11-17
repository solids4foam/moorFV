/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "sixDoFRigidBodyMotionAxisFvBeamConstraint.H"
#include "addToRunTimeSelectionTable.H"
#include "sixDoFRigidBodyMotionFvBeam.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sixDoFRigidBodyMotionFvBeamConstraints
{
    defineTypeNameAndDebug(axis, 0);

    addToRunTimeSelectionTable
    (
        sixDoFRigidBodyMotionFvBeamConstraint,
        axis,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionFvBeamConstraints::axis::axis
(
    const word& name,
    const dictionary& sDoFRBMCDict,
    const sixDoFRigidBodyMotionFvBeam& motion
)
:
    sixDoFRigidBodyMotionFvBeamConstraint(name, sDoFRBMCDict, motion),
    axis_()
{
    read(sDoFRBMCDict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionFvBeamConstraints::axis::~axis()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::sixDoFRigidBodyMotionFvBeamConstraints::axis::constrainTranslation
(
    pointConstraint& pc
) const
{}


void Foam::sixDoFRigidBodyMotionFvBeamConstraints::axis::constrainRotation
(
    pointConstraint& pc
) const
{
    pc.combine(pointConstraint(Tuple2<label, vector>(2, axis_)));
}


bool Foam::sixDoFRigidBodyMotionFvBeamConstraints::axis::read
(
    const dictionary& sDoFRBMCDict
)
{
    sixDoFRigidBodyMotionFvBeamConstraint::read(sDoFRBMCDict);

    sDoFRBMCCoeffs_.readEntry("axis", axis_);

    scalar magFixedAxis(mag(axis_));

    if (magFixedAxis > VSMALL)
    {
        axis_ /= magFixedAxis;
    }
    else
    {
        FatalErrorInFunction
            << "axis has zero length"
            << abort(FatalError);
    }

    return true;
}


void Foam::sixDoFRigidBodyMotionFvBeamConstraints::axis::write
(
    Ostream& os
) const
{
    os.writeEntry("axis", axis_);
}

// ************************************************************************* //
