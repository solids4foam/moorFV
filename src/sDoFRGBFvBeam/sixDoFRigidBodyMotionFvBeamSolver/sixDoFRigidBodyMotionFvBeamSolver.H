/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2017 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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
    Foam::sixDoFRigidBodyMotionFvBeamSolver

Group
   grpMeshMotionSolvers

Description
    6-DoF solid-body mesh motion solver for an fvMesh.

    Applies SLERP interpolation of movement as function of distance to
    the object surface.

SourceFiles
    sixDoFRigidBodyMotionFvBeamSolver.C

\*---------------------------------------------------------------------------*/

#ifndef sixDoFRigidBodyMotionFvBeamSolver_H
#define sixDoFRigidBodyMotionFvBeamSolver_H

#include "displacementMotionSolver.H"
#include "sixDoFRigidBodyMotionFvBeam.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class sixDoFRigidBodyMotionFvBeamSolver Declaration
\*---------------------------------------------------------------------------*/

class sixDoFRigidBodyMotionFvBeamSolver
:
    public displacementMotionSolver
{
    // Private data

        //- Six DoF motion object
        sixDoFRigidBodyMotionFvBeam motion_;

        wordRes patches_;

        //- Patches to integrate forces
        const labelHashSet patchSet_;

        //- Inner morphing distance (limit of solid-body region)
        const scalar di_;

        //- Outer morphing distance (limit of linear interpolation region)
        const scalar do_;

        //- Test-mode in which only the gravitational body-force is applied
        bool test_;

        //- Reference density required by the forces object for
        //  incompressible calculations, required if rho == rhoInf
        scalar rhoInf_;

        //- Name of density field, optional unless used for an
        //  incompressible simulation, when this needs to be specified
        //  as rhoInf
        word rhoName_;

        //- Current interpolation scale (1 at patches, 0 at distance_)
        pointScalarField scale_;

        //- Current time index (used for updating)
        label curTimeIndex_;

        //- Name of the uniformVectorField for CofG velocity
        word CofGvelocity_;


    // Private Member Functions

        //- No copy construct
        sixDoFRigidBodyMotionFvBeamSolver
        (
            const sixDoFRigidBodyMotionFvBeamSolver&
        ) = delete;

        //- No copy assignment
        void operator=(const sixDoFRigidBodyMotionFvBeamSolver&) = delete;


public:

    //- Runtime type information
    TypeName("sixDoFRigidBodyMotionFvBeam");


    // Constructors

        //- Construct from polyMesh and IOdictionary
        sixDoFRigidBodyMotionFvBeamSolver
        (
            const polyMesh&,
            const IOdictionary& dict
        );


    //- Destructor
    ~sixDoFRigidBodyMotionFvBeamSolver() = default;


    // Member Functions

        //- Return the six DoF motion object
        const sixDoFRigidBodyMotionFvBeam& motion() const;

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> curPoints() const;

        //- Solve for motion
        virtual void solve();

        //- Write state using stream options
        virtual bool writeObject
        (
            IOstreamOption streamOpt,
            const bool valid
        ) const;

        //- Read dynamicMeshDict dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
