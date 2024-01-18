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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::point&
Foam::sixDoFRigidBodyMotionFvBeamState::centreOfRotation() const
{
    return centreOfRotation_;
}


inline const Foam::tensor& Foam::sixDoFRigidBodyMotionFvBeamState::Q() const
{
    return Q_;
}


inline const Foam::vector& Foam::sixDoFRigidBodyMotionFvBeamState::v() const
{
    return v_;
}


inline const Foam::vector& Foam::sixDoFRigidBodyMotionFvBeamState::a() const
{
    return a_;
}


inline const Foam::vector& Foam::sixDoFRigidBodyMotionFvBeamState::pi() const
{
    return pi_;
}


inline const Foam::vector& Foam::sixDoFRigidBodyMotionFvBeamState::tau() const
{
    return tau_;
}


inline Foam::point& Foam::sixDoFRigidBodyMotionFvBeamState::centreOfRotation()
{
    return centreOfRotation_;
}


inline Foam::tensor& Foam::sixDoFRigidBodyMotionFvBeamState::Q()
{
    return Q_;
}


inline Foam::vector& Foam::sixDoFRigidBodyMotionFvBeamState::v()
{
    return v_;
}


inline Foam::vector& Foam::sixDoFRigidBodyMotionFvBeamState::a()
{
    return a_;
}


inline Foam::vector& Foam::sixDoFRigidBodyMotionFvBeamState::pi()
{
    return pi_;
}


inline Foam::vector& Foam::sixDoFRigidBodyMotionFvBeamState::tau()
{
    return tau_;
}


// ************************************************************************* //
