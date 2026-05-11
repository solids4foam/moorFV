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

#include "sixDoFRigidBodyBeamMotionOrientationConstraint.H"
#include "addToRunTimeSelectionTable.H"
#include "sixDoFRigidBodyBeamMotion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sixDoFRigidBodyBeamMotionConstraints
{
    defineTypeNameAndDebug(orientation, 0);

    addToRunTimeSelectionTable
    (
        sixDoFRigidBodyBeamMotionConstraint,
        orientation,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyBeamMotionConstraints::orientation::orientation
(
    const word& name,
    const dictionary& sDoFRBMCDict,
    const sixDoFRigidBodyBeamMotion& motion
)
:
    sixDoFRigidBodyBeamMotionConstraint(name, sDoFRBMCDict, motion)
{
    read(sDoFRBMCDict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyBeamMotionConstraints::orientation::~orientation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::sixDoFRigidBodyBeamMotionConstraints::orientation::
constrainTranslation
(
    pointConstraint& pc
) const
{}


void Foam::sixDoFRigidBodyBeamMotionConstraints::orientation::
constrainRotation
(
    pointConstraint& pc
) const
{
    pc.combine(pointConstraint(Tuple2<label, vector>(3, Zero)));
}


bool Foam::sixDoFRigidBodyBeamMotionConstraints::orientation::read
(
    const dictionary& sDoFRBMCDict
)
{
    sixDoFRigidBodyBeamMotionConstraint::read(sDoFRBMCDict);

    return true;
}


void Foam::sixDoFRigidBodyBeamMotionConstraints::orientation::write
(
    Ostream& os
) const
{
}

// ************************************************************************* //
