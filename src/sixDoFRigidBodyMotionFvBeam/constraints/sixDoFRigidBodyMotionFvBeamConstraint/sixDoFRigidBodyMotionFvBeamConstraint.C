/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "sixDoFRigidBodyMotionFvBeamConstraint.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sixDoFRigidBodyMotionFvBeamConstraint, 0);
    defineRunTimeSelectionTable(sixDoFRigidBodyMotionFvBeamConstraint, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionFvBeamConstraint::sixDoFRigidBodyMotionFvBeamConstraint
(
    const word& name,
    const dictionary& sDoFRBMCDict,
    const sixDoFRigidBodyMotionFvBeam& motion
)
:
    name_(name),
    sDoFRBMCCoeffs_(sDoFRBMCDict),
    motion_(motion)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionFvBeamConstraint::~sixDoFRigidBodyMotionFvBeamConstraint()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::sixDoFRigidBodyMotionFvBeamConstraint::read
(
    const dictionary& sDoFRBMCDict
)
{
    sDoFRBMCCoeffs_ = sDoFRBMCDict;

    return true;
}


void Foam::sixDoFRigidBodyMotionFvBeamConstraint::write(Ostream& os) const
{}


// ************************************************************************* //
