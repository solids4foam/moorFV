/*---------------------------------------------------------------------------* \
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "sixDoFRigidBodyBeamMotionRestraint.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sixDoFRigidBodyBeamMotionRestraint, 0);
    defineRunTimeSelectionTable(sixDoFRigidBodyBeamMotionRestraint, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyBeamMotionRestraint::sixDoFRigidBodyBeamMotionRestraint
(
    const word& name,
    const dictionary& sDoFRBMRDict,
    const Time& time
)
:
    name_(name),
    sDoFRBMRCoeffs_(sDoFRBMRDict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyBeamMotionRestraint::~sixDoFRigidBodyBeamMotionRestraint()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::sixDoFRigidBodyBeamMotionRestraint::read
(
    const dictionary& sDoFRBMRDict
)
{
    sDoFRBMRCoeffs_ = sDoFRBMRDict;

    return true;
}


// ************************************************************************* //
