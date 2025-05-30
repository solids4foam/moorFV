/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "sixDoFRigidBodyMotionFvBeamState.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionFvBeamState::sixDoFRigidBodyMotionFvBeamState()
:
    centreOfRotation_(Zero),
    Q_(I),
    v_(Zero),
    a_(Zero),
    pi_(Zero),
    tau_(Zero)
{}


Foam::sixDoFRigidBodyMotionFvBeamState::sixDoFRigidBodyMotionFvBeamState
(
    const dictionary& dict
)
:
    centreOfRotation_
    (
        dict.getOrDefault
        (
            "centreOfRotation",
            dict.getOrDefault("centreOfMass", vector::zero)
        )
    ),
    Q_(dict.getOrDefault("orientation", tensor::I)),
    v_(dict.getOrDefault("velocity", vector::zero)),
    a_(dict.getOrDefault("acceleration", vector::zero)),
    pi_(dict.getOrDefault("angularMomentum", vector::zero)),
    tau_(dict.getOrDefault("torque", vector::zero))
{}


// ************************************************************************* //
