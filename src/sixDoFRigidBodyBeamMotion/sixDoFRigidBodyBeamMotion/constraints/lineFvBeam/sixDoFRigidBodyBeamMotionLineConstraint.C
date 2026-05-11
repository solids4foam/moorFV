/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "sixDoFRigidBodyBeamMotionLineConstraint.H"
#include "addToRunTimeSelectionTable.H"
#include "sixDoFRigidBodyBeamMotion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sixDoFRigidBodyBeamMotionConstraints
{
    defineTypeNameAndDebug(line, 0);

    addToRunTimeSelectionTable
    (
        sixDoFRigidBodyBeamMotionConstraint,
        line,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyBeamMotionConstraints::line::line
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

Foam::sixDoFRigidBodyBeamMotionConstraints::line::~line()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::sixDoFRigidBodyBeamMotionConstraints::line::setCentreOfRotation
(
    point& CofR
) const
{
    CofR = centreOfRotation_;
}


void Foam::sixDoFRigidBodyBeamMotionConstraints::line::constrainTranslation
(
    pointConstraint& pc
) const
{
    pc.combine(pointConstraint(Tuple2<label, vector>(2, direction_)));
}


void Foam::sixDoFRigidBodyBeamMotionConstraints::line::constrainRotation
(
    pointConstraint& pc
) const
{}


bool Foam::sixDoFRigidBodyBeamMotionConstraints::line::read
(
    const dictionary& sDoFRBMCDict
)
{
    sixDoFRigidBodyBeamMotionConstraint::read(sDoFRBMCDict);

    centreOfRotation_ = sDoFRBMCCoeffs_.getOrDefault
    (
        "centreOfRotation",
        motion_.initialCentreOfMass()
    );

    sDoFRBMCCoeffs_.readEntry("direction", direction_);

    scalar magDir(mag(direction_));

    if (magDir > VSMALL)
    {
        direction_ /= magDir;
    }
    else
    {
        FatalErrorInFunction
            << "line direction has zero length"
            << abort(FatalError);
    }

    return true;
}


void Foam::sixDoFRigidBodyBeamMotionConstraints::line::write
(
    Ostream& os
) const
{
    os.writeEntry("centreOfRotation", centreOfRotation_);
    os.writeEntry("direction", direction_);
}

// ************************************************************************* //
