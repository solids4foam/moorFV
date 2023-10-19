/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2016-2019 OpenCFD Ltd.
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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline Foam::tensor Foam::sixDoFRigidBodyMotionFvBeam::rotationTensorX
(
    scalar phi
) const
{
    return tensor
    (
        1, 0, 0,
        0, Foam::cos(phi), -Foam::sin(phi),
        0, Foam::sin(phi), Foam::cos(phi)
    );
}


inline Foam::tensor Foam::sixDoFRigidBodyMotionFvBeam::rotationTensorY
(
    scalar phi
) const
{
    return tensor
    (
        Foam::cos(phi), 0, Foam::sin(phi),
        0, 1, 0,
        -Foam::sin(phi), 0, Foam::cos(phi)
    );
}


inline Foam::tensor Foam::sixDoFRigidBodyMotionFvBeam::rotationTensorZ
(
    scalar phi
) const
{
    return tensor
    (
        Foam::cos(phi), -Foam::sin(phi), 0,
        Foam::sin(phi), Foam::cos(phi), 0,
        0, 0, 1
    );
}


inline Foam::Tuple2<Foam::tensor, Foam::vector>
Foam::sixDoFRigidBodyMotionFvBeam::rotate
(
    const tensor& Q0,
    const vector& pi0,
    const scalar deltaT
) const
{
    Tuple2<tensor, vector> Qpi(Q0, pi0);
    tensor& Q = Qpi.first();
    vector& pi = Qpi.second();

    tensor R = rotationTensorX(0.5*deltaT*pi.x()/momentOfInertia_.xx());
    pi = pi & R;
    Q = Q & R;

    R = rotationTensorY(0.5*deltaT*pi.y()/momentOfInertia_.yy());
    pi = pi & R;
    Q = Q & R;

    R = rotationTensorZ(deltaT*pi.z()/momentOfInertia_.zz());
    pi = pi & R;
    Q = Q & R;

    R = rotationTensorY(0.5*deltaT*pi.y()/momentOfInertia_.yy());
    pi = pi & R;
    Q = Q & R;

    R = rotationTensorX(0.5*deltaT*pi.x()/momentOfInertia_.xx());
    pi = pi & R;
    Q = Q & R;

    return Qpi;
}


inline const Foam::PtrList<Foam::sixDoFRigidBodyMotionFvBeamRestraint>&
Foam::sixDoFRigidBodyMotionFvBeam::restraints() const
{
    return restraints_;
}


inline const Foam::PtrList<Foam::sixDoFRigidBodyMotionFvBeamConstraint>&
Foam::sixDoFRigidBodyMotionFvBeam::constraints() const
{
    return constraints_;
}


inline const Foam::point&
Foam::sixDoFRigidBodyMotionFvBeam::initialCentreOfRotation() const
{
    return initialCentreOfRotation_;
}


inline const Foam::tensor&
Foam::sixDoFRigidBodyMotionFvBeam::initialQ() const
{
    return initialQ_;
}


inline const Foam::tensor& Foam::sixDoFRigidBodyMotionFvBeam::Q() const
{
    return motionState_.Q();
}


inline const Foam::vector& Foam::sixDoFRigidBodyMotionFvBeam::v() const
{
    return motionState_.v();
}


inline const Foam::vector& Foam::sixDoFRigidBodyMotionFvBeam::a() const
{
    return motionState_.a();
}


inline const Foam::vector& Foam::sixDoFRigidBodyMotionFvBeam::pi() const
{
    return motionState_.pi();
}


inline const Foam::vector& Foam::sixDoFRigidBodyMotionFvBeam::tau() const
{
    return motionState_.tau();
}


inline Foam::point& Foam::sixDoFRigidBodyMotionFvBeam::initialCentreOfRotation()
{
    return initialCentreOfRotation_;
}


inline Foam::tensor& Foam::sixDoFRigidBodyMotionFvBeam::initialQ()
{
    return initialQ_;
}


inline Foam::tensor& Foam::sixDoFRigidBodyMotionFvBeam::Q()
{
    return motionState_.Q();
}


inline Foam::vector& Foam::sixDoFRigidBodyMotionFvBeam::v()
{
    return motionState_.v();
}


inline Foam::vector& Foam::sixDoFRigidBodyMotionFvBeam::a()
{
    return motionState_.a();
}


inline Foam::vector& Foam::sixDoFRigidBodyMotionFvBeam::pi()
{
    return motionState_.pi();
}


inline Foam::vector& Foam::sixDoFRigidBodyMotionFvBeam::tau()
{
    return motionState_.tau();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

inline Foam::scalar Foam::sixDoFRigidBodyMotionFvBeam::mass() const
{
    return mass_;
}


inline const Foam::diagTensor&
Foam::sixDoFRigidBodyMotionFvBeam::momentOfInertia() const
{
    return momentOfInertia_;
}


inline const Foam::sixDoFRigidBodyMotionFvBeamState&
Foam::sixDoFRigidBodyMotionFvBeam::state() const
{
    return motionState_;
}


inline const Foam::point& Foam::sixDoFRigidBodyMotionFvBeam::centreOfRotation() const
{
    return motionState_.centreOfRotation();
}


inline const Foam::point&
Foam::sixDoFRigidBodyMotionFvBeam::initialCentreOfMass() const
{
    return initialCentreOfMass_;
}


inline Foam::point Foam::sixDoFRigidBodyMotionFvBeam::centreOfMass() const
{
    return transform(initialCentreOfMass_);
}


inline Foam::vector Foam::sixDoFRigidBodyMotionFvBeam::momentArm() const
{
    return centreOfMass() - motionState_.centreOfRotation();
}


inline const Foam::tensor&
Foam::sixDoFRigidBodyMotionFvBeam::orientation() const
{
    return Q();
}


inline Foam::vector Foam::sixDoFRigidBodyMotionFvBeam::omega() const
{
    return  Q() & (inv(momentOfInertia_) & pi());
}

inline const Foam::Time& Foam::sixDoFRigidBodyMotionFvBeam::time() const
{
    return time_;
}

inline bool Foam::sixDoFRigidBodyMotionFvBeam::report() const
{
    return report_;
}


inline void Foam::sixDoFRigidBodyMotionFvBeam::newTime()
{
    motionState0_ = motionState_;
}


inline Foam::point& Foam::sixDoFRigidBodyMotionFvBeam::centreOfRotation()
{
    return motionState_.centreOfRotation();
}


inline Foam::point Foam::sixDoFRigidBodyMotionFvBeam::velocity
(
    const point& pt
) const
{
    return (omega() ^ (pt - centreOfRotation())) + v();
}


inline Foam::point Foam::sixDoFRigidBodyMotionFvBeam::transform
(
    const point& initialPoint
) const
{
    return
    (
        centreOfRotation()
      + (Q() & initialQ_.T() & (initialPoint - initialCentreOfRotation_))
    );
}

// ************************************************************************* //
