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

Class
    Foam::sixDoFRigidBodyMotionConstraints::plane

Description
    Translation constraint on the centre of rotation:
        may only move on a plane.

    If 'centreOfRotation' is not provided in the dictionary the centre of mass
    is used.

SourceFiles
    sixDoFRigidBodyMotionPlaneConstraint.C

\*---------------------------------------------------------------------------*/

#ifndef sixDoFRigidBodyMotionPlaneFvBeamConstraint_H
#define sixDoFRigidBodyMotionPlaneFvBeamConstraint_H

#include "sixDoFRigidBodyMotionFvBeamConstraint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace sixDoFRigidBodyMotionFvBeamConstraints
{

/*---------------------------------------------------------------------------*\
                          Class plane Declaration
\*---------------------------------------------------------------------------*/

class plane
:
    public sixDoFRigidBodyMotionFvBeamConstraint
{
    // Private data

        //- Centre of rotation on plane
        point centreOfRotation_;

        //- Normal to plane
        vector normal_;


public:

    //- Runtime type information
    TypeName("plane");


    // Constructors

        //- Construct from components
        plane
        (
            const word& name,
            const dictionary& sDoFRBMCDict,
            const sixDoFRigidBodyMotionFvBeam& motion
        );

        //- Construct and return a clone
        virtual autoPtr<sixDoFRigidBodyMotionFvBeamConstraint> clone() const
        {
            return autoPtr<sixDoFRigidBodyMotionFvBeamConstraint>
            (
                new plane(*this)
            );
        }


    //- Destructor
    virtual ~plane();


    // Member Functions

        //- Set the centre of rotation to the projection of the
        //  centre of mass onto the plane
        virtual void setCentreOfRotation(point&) const;

        //- Apply and accumulate translational constraints
        virtual void constrainTranslation(pointConstraint&) const;

        //- Apply and accumulate rotational constraints
        virtual void constrainRotation(pointConstraint&) const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& sDoFRBMCCoeff);

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidBodyMotionFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
