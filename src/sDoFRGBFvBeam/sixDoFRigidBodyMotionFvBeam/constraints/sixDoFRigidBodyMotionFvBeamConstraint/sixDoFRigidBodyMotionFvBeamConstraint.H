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

Namespace
    Foam::sixDoFRigidBodyMotionFvBeamConstraints

Description
    Namespace for six DoF motion constraints


Class
    Foam::sixDoFRigidBodyMotionFvBeamConstraint

Description
    Base class for defining constraints for sixDoF motions

SourceFiles
    sixDoFRigidBodyMotionFvBeamConstraint.C
    sixDoFRigidBodyMotionFvBeamConstraintNew.C

\*---------------------------------------------------------------------------*/

#ifndef sixDoFRigidBodyMotionFvBeamConstraint_H
#define sixDoFRigidBodyMotionFvBeamConstraint_H

#include "Time.H"
#include "dictionary.H"
#include "autoPtr.H"
#include "point.H"
#include "pointConstraint.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class sixDoFRigidBodyMotionFvBeam;


/*---------------------------------------------------------------------------*\
                Class sixDoFRigidBodyMotionFvBeamConstraint Declaration
\*---------------------------------------------------------------------------*/

class sixDoFRigidBodyMotionFvBeamConstraint
{

protected:

    // Protected data

        //- Name of the constraint
        word name_;

        //- Constraint model specific coefficient dictionary
        dictionary sDoFRBMCCoeffs_;

        //- Reference to the body motion
        const sixDoFRigidBodyMotionFvBeam& motion_;


public:

    //- Runtime type information
    TypeName("sixDoFRigidBodyMotionFvBeamConstraint");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            sixDoFRigidBodyMotionFvBeamConstraint,
            dictionary,
            (
                const word& name,
                const dictionary& sDoFRBMCDict,
                const sixDoFRigidBodyMotionFvBeam& motion
            ),
            (name, sDoFRBMCDict, motion)
        );


    // Constructors

        //- Construct from the sDoFRBMCDict dictionary and Time
        sixDoFRigidBodyMotionFvBeamConstraint
        (
            const word& name,
            const dictionary& sDoFRBMCDict,
            const sixDoFRigidBodyMotionFvBeam& motion
        );

        //- Construct and return a clone
        virtual autoPtr<sixDoFRigidBodyMotionFvBeamConstraint> clone() const = 0;


    // Selectors

        //- Select constructed from the sDoFRBMCDict dictionary and Time
        static autoPtr<sixDoFRigidBodyMotionFvBeamConstraint> New
        (
            const word& name,
            const dictionary& sDoFRBMCDict,
            const sixDoFRigidBodyMotionFvBeam& motion
        );


    //- Destructor
    virtual ~sixDoFRigidBodyMotionFvBeamConstraint();


    // Member Functions

        //- Return the name
        const word& name() const
        {
            return name_;
        }

        //- Set the centre of rotation if not the centre of mass
        virtual void setCentreOfRotation(point&) const
        {}

        //- Apply and accumulate translational constraints
        virtual void constrainTranslation(pointConstraint&) const = 0;

        //- Apply and accumulate rotational constraints
        virtual void constrainRotation(pointConstraint&) const = 0;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& sDoFRBMCDict);

        // Access

            // Return access to sDoFRBMCCoeffs
            inline const dictionary& coeffDict() const
            {
                return sDoFRBMCCoeffs_;
            }

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
