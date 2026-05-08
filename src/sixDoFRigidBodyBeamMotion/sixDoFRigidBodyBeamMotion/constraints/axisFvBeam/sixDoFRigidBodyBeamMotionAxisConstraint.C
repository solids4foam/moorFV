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

#include "sixDoFRigidBodyBeamMotionAxisConstraint.H"
#include "addToRunTimeSelectionTable.H"
#include "sixDoFRigidBodyBeamMotion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sixDoFRigidBodyBeamMotionConstraints
{
    defineTypeNameAndDebug(axis, 0);

    addToRunTimeSelectionTable
    (
        sixDoFRigidBodyBeamMotionConstraint,
        axis,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyBeamMotionConstraints::axis::axis
(
    const word& name,
    const dictionary& sDoFRBMCDict,
    const sixDoFRigidBodyBeamMotion& motion
)
:
    sixDoFRigidBodyBeamMotionConstraint(name, sDoFRBMCDict, motion),
    axis_()
{
    read(sDoFRBMCDict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyBeamMotionConstraints::axis::~axis()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::sixDoFRigidBodyBeamMotionConstraints::axis::constrainTranslation
(
    pointConstraint& pc
) const
{}


void Foam::sixDoFRigidBodyBeamMotionConstraints::axis::constrainRotation
(
    pointConstraint& pc
) const
{
    pc.combine(pointConstraint(Tuple2<label, vector>(2, axis_)));
}


bool Foam::sixDoFRigidBodyBeamMotionConstraints::axis::read
(
    const dictionary& sDoFRBMCDict
)
{
    sixDoFRigidBodyBeamMotionConstraint::read(sDoFRBMCDict);

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


void Foam::sixDoFRigidBodyBeamMotionConstraints::axis::write
(
    Ostream& os
) const
{
    os.writeEntry("axis", axis_);
}

// ************************************************************************* //
