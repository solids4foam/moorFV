/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
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
    Foam::sixDoFFvBeamSolver

Group
    grpSixDoFRigidBodySolvers

Description

SourceFiles
    sixDoFFvBeamSolver.C
    newsixDoFFvBeamSolver.C

\*---------------------------------------------------------------------------*/

#ifndef sixDoFFvBeamSolver_H
#define sixDoFFvBeamSolver_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "sixDoFRigidBodyMotionFvBeam.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class sixDoFFvBeamSolver Declaration
\*---------------------------------------------------------------------------*/

class sixDoFFvBeamSolver
{
protected:

    // Protected data

        //- The rigid body
        sixDoFRigidBodyMotionFvBeam& body_;

        //- Model dictionary
        dictionary dict_;


    // Protected member functions

        //- Return the current centre of rotation
        inline point& centreOfRotation();

        //- Return the orientation
        inline tensor& Q();

        //- Return non-const access to vector
        inline vector& v();

        //- Return non-const access to acceleration
        inline vector& a();

        //- Return non-const access to angular momentum
        inline vector& pi();

        //- Return non-const access to torque
        inline vector& tau();

        //- Return the centre of rotation at previous time-step
        inline const point& centreOfRotation0() const;

        //- Return the orientation at previous time-step
        inline const tensor& Q0() const;

        //- Return the velocity at previous time-step
        inline const vector& v0() const;

        //- Return the acceleration at previous time-step
        inline const vector& a0() const;

        //- Return the angular momentum at previous time-step
        inline const vector& pi0() const;

        //- Return the torque at previous time-step
        inline const vector& tau0() const;

        //- Acceleration damping coefficient (for steady-state simulations)
        inline scalar aDamp() const;

        //- Translational constraint tensor
        inline tensor tConstraints() const;

        //- Rotational constraint tensor
        inline tensor rConstraints() const;

        //- Apply rotation tensors to Q0 for the given torque (pi) and deltaT
        //  and return the rotated Q and pi as a tuple
        inline Tuple2<tensor, vector> rotate
        (
            const tensor& Q0,
            const vector& pi,
            const scalar deltaT
        ) const;

        //- Update and relax accelerations from the force and torque
        inline void updateAcceleration
        (
            const vector& fGlobal,
            const vector& tauGlobal
        );


public:

    //- Runtime type information
    TypeName("sixDoFFvBeamSolver");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            sixDoFFvBeamSolver,
            dictionary,
            (
                const dictionary& dict,
                sixDoFRigidBodyMotionFvBeam& body
            ),
            (dict, body)
        );


    // Constructors

        // Construct for given body
        sixDoFFvBeamSolver(const dictionary& dict, sixDoFRigidBodyMotionFvBeam& body);

        //- Construct and return a clone
        virtual autoPtr<sixDoFFvBeamSolver> clone() const = 0;


    //- Destructor
    virtual ~sixDoFFvBeamSolver();


    // Selectors

        static autoPtr<sixDoFFvBeamSolver> New
        (
            const dictionary& dict,
            sixDoFRigidBodyMotionFvBeam& body
        );


    // Member Functions

        //- Drag coefficient
        virtual void solve
        (
            bool firstIter,
            const vector& fGlobal,
            const vector& tauGlobal,
            scalar deltaT,
            scalar deltaT0
        ) = 0;


    //- Write
    void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "sixDoFFvBeamSolverI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
