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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline Foam::tensor
Foam::sixDoFRigidBodyMotionNew::rotationTensorX(scalar phi) const
{
    return tensor
    (
        1, 0, 0,
        0, Foam::cos(phi), -Foam::sin(phi),
        0, Foam::sin(phi), Foam::cos(phi)
    );
}


inline Foam::tensor
Foam::sixDoFRigidBodyMotionNew::rotationTensorY(scalar phi) const
{
    return tensor
    (
        Foam::cos(phi), 0, Foam::sin(phi),
        0, 1, 0,
        -Foam::sin(phi), 0, Foam::cos(phi)
    );
}


inline Foam::tensor
Foam::sixDoFRigidBodyMotionNew::rotationTensorZ(scalar phi) const
{
    return tensor
    (
        Foam::cos(phi), -Foam::sin(phi), 0,
        Foam::sin(phi), Foam::cos(phi), 0,
        0, 0, 1
    );
}


inline void Foam::sixDoFRigidBodyMotionNew::rotate
(
    tensor& Q,
    vector& pi,
    scalar deltaT
) const
{
    tensor R;

    R = rotationTensorX(0.5*deltaT*pi.x()/momentOfInertia_.xx());
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
}


inline const Foam::sixDoFRigidBodyMotionState&
Foam::sixDoFRigidBodyMotionNew::motionState() const
{
    return motionState_;
}


inline const Foam::PtrList<Foam::sixDoFRigidBodyMotionRestraintNew>&
Foam::sixDoFRigidBodyMotionNew::restraints() const
{
    return restraints_;
}


inline const Foam::wordList& Foam::sixDoFRigidBodyMotionNew::restraintNames() const
{
    return restraintNames_;
}


inline const Foam::PtrList<Foam::sixDoFRigidBodyMotionConstraintNew>&
Foam::sixDoFRigidBodyMotionNew::constraints() const
{
    return constraints_;
}


inline const Foam::wordList&
Foam::sixDoFRigidBodyMotionNew::constraintNames() const
{
    return constraintNames_;
}


inline Foam::label Foam::sixDoFRigidBodyMotionNew::maxConstraintIterations() const
{
    return maxConstraintIterations_;
}


inline const Foam::point&
Foam::sixDoFRigidBodyMotionNew::initialCentreOfMass() const
{
    return initialCentreOfMass_;
}


inline const Foam::tensor&
Foam::sixDoFRigidBodyMotionNew::initialQ() const
{
    return initialQ_;
}


inline const Foam::tensor& Foam::sixDoFRigidBodyMotionNew::Q() const
{
    return motionState_.Q();
}


inline const Foam::vector& Foam::sixDoFRigidBodyMotionNew::v() const
{
    return motionState_.v();
}


inline const Foam::vector& Foam::sixDoFRigidBodyMotionNew::a() const
{
    return motionState_.a();
}


inline const Foam::vector& Foam::sixDoFRigidBodyMotionNew::pi() const
{
    return motionState_.pi();
}


inline const Foam::vector& Foam::sixDoFRigidBodyMotionNew::tau() const
{
    return motionState_.tau();
}


inline Foam::point& Foam::sixDoFRigidBodyMotionNew::initialCentreOfMass()
{
    return initialCentreOfMass_;
}


inline Foam::tensor& Foam::sixDoFRigidBodyMotionNew::initialQ()
{
    return initialQ_;
}


inline Foam::tensor& Foam::sixDoFRigidBodyMotionNew::Q()
{
    return motionState_.Q();
}


inline Foam::vector& Foam::sixDoFRigidBodyMotionNew::v()
{
    return motionState_.v();
}


inline Foam::vector& Foam::sixDoFRigidBodyMotionNew::a()
{
    return motionState_.a();
}


inline Foam::vector& Foam::sixDoFRigidBodyMotionNew::pi()
{
    return motionState_.pi();
}


inline Foam::vector& Foam::sixDoFRigidBodyMotionNew::tau()
{
    return motionState_.tau();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

inline Foam::tmp<Foam::pointField>
Foam::sixDoFRigidBodyMotionNew::currentPosition(const pointField& pInitial) const
{
    return
    (
        centreOfMass()
      + (Q() & initialQ_.T() & (pInitial - initialCentreOfMass_)))
    ;
}


inline Foam::point Foam::sixDoFRigidBodyMotionNew::currentPosition
(
    const point& pInitial
) const
{
    return
    (
        centreOfMass()
      + (Q() & initialQ_.T() & (pInitial - initialCentreOfMass_))
    );
}


inline Foam::vector Foam::sixDoFRigidBodyMotionNew::currentOrientation
(
    const vector& vInitial
) const
{
    return (Q() & initialQ_.T() & vInitial);
}


inline const Foam::tensor&
Foam::sixDoFRigidBodyMotionNew::orientation() const
{
    return Q();
}


inline Foam::vector Foam::sixDoFRigidBodyMotionNew::omega() const
{
    return  Q() & (inv(momentOfInertia_) & pi());
}


inline Foam::point Foam::sixDoFRigidBodyMotionNew::currentVelocity
(
    const point& pt
) const
{
    return (omega() ^ (pt - centreOfMass())) + v();
}


inline const Foam::point& Foam::sixDoFRigidBodyMotionNew::centreOfMass() const
{
    return motionState_.centreOfMass();
}


inline const Foam::diagTensor&
Foam::sixDoFRigidBodyMotionNew::momentOfInertia() const
{
    return momentOfInertia_;
}


inline Foam::scalar Foam::sixDoFRigidBodyMotionNew::mass() const
{
    return mass_;
}


inline bool Foam::sixDoFRigidBodyMotionNew::report() const
{
    return report_;
}


inline Foam::point& Foam::sixDoFRigidBodyMotionNew::centreOfMass()
{
    return motionState_.centreOfMass();
}


inline Foam::diagTensor& Foam::sixDoFRigidBodyMotionNew::momentOfInertia()
{
    return momentOfInertia_;
}


inline Foam::scalar& Foam::sixDoFRigidBodyMotionNew::mass()
{
    return mass_;
}

// ************************************************************************* //
