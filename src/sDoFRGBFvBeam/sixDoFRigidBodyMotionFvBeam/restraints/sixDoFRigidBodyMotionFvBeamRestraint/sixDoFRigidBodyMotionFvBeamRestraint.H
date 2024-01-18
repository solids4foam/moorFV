/*---------------------------------------------------------------------------*\
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

Namespace
    Foam::sixDoFRigidBodyMotionFvBeamRestraints

Description
    Namespace for six DoF motion restraints


Class
    Foam::sixDoFRigidBodyMotionFvBeamRestraint

Description
    Base class for defining restraints for sixDoF motions

SourceFiles
    sixDoFRigidBodyMotionFvBeamRestraint.C
    sixDoFRigidBodyMotionFvBeamRestraintNew.C

\*---------------------------------------------------------------------------*/

#ifndef sixDoFRigidBodyMotionFvBeamRestraint_H
#define sixDoFRigidBodyMotionFvBeamRestraint_H

#include "Time.H"
#include "dictionary.H"
#include "autoPtr.H"
#include "vector.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class sixDoFRigidBodyMotionFvBeam;


/*---------------------------------------------------------------------------*\
                Class sixDoFRigidBodyMotionFvBeamRestraint Declaration
\*---------------------------------------------------------------------------*/

class sixDoFRigidBodyMotionFvBeamRestraint
{

protected:

    // Protected data

        //- Name of the restraint
        word name_;

        //- Restraint model specific coefficient dictionary
        dictionary sDoFRBMRCoeffs_;


public:

    //- Runtime type information
    TypeName("sixDoFRigidBodyMotionFvBeamRestraint");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            sixDoFRigidBodyMotionFvBeamRestraint,
            dictionary,
            (const word& name, const dictionary& sDoFRBMRDict, const Time& time),
            (name, sDoFRBMRDict, time)
        );


    // Constructors

        //- Construct from the sDoFRBMRDict dictionary and Time
        sixDoFRigidBodyMotionFvBeamRestraint
        (
            const word& name,
            const dictionary& sDoFRBMRDict,
            const Time& time
        );

        //- Construct and return a clone
        virtual autoPtr<sixDoFRigidBodyMotionFvBeamRestraint> clone() const = 0;


    // Selectors

        //- Select constructed from the sDoFRBMRDict dictionary and Time
        static autoPtr<sixDoFRigidBodyMotionFvBeamRestraint> New
        (
            const word& name,
            const dictionary& sDoFRBMRDict,
            const Time& time
        );


    //- Destructor
    virtual ~sixDoFRigidBodyMotionFvBeamRestraint();


    // Member Functions

        //- Return the name
        const word& name() const
        {
            return name_;
        }

        //- Calculate the restraint position, force and moment.
        //  Global reference frame vectors.
        virtual void restrain
        (
            const sixDoFRigidBodyMotionFvBeam& motion,
            vector& restraintPosition,
            vector& restraintForce,
            vector& restraintMoment
        ) const = 0;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& sDoFRBMRDict);

        // Access

            // Return access to sDoFRBMRCoeffs
            inline const dictionary& coeffDict() const
            {
                return sDoFRBMRCoeffs_;
            }

        //- Write
        virtual void write(Ostream&) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
