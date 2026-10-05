/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
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
    Foam::sixDoFFvBeamSolvers::beamFoamCoupled

Group
    grpSixDoFRigidBodySolvers

Description
    Rigid body solved monolithically with its finiteVolumeBeam mooring line
    in beamFoam's BlockEigen system.

    The fluid force and moment (with the weight) and any other restraints,
    evaluated at the current state, are passed to the finiteVolumeBeam
    restraint as an applied load. beamFoam solves the beam and the body's
    translation and rotation in one Newton iteration, and the accepted body
    state is copied into the motion as it is: there is no separate Newmark
    update and no acceleration relaxation.

    Motion constraints are applied inside the monolithic solve: the body's
    Newton increments are projected onto the free directions.

    Requires one finiteVolumeBeam restraint and the centre of rotation at the
    centre of mass.

Usage
    \verbatim
    solver
    {
        type            beamFoamCoupled;
        beta            0.25;   // Newmark beta (default 0.25)
        gamma           0.5;    // Newmark gamma (default 0.5)
        linearDamping   0;      // viscous force -c v [N s/m]
        angularDamping  0;      // viscous moment -c omega [N m s]
        jacobianCheck   0;      // finite-difference check, first N steps
    }
    \endverbatim

SourceFiles
    beamFoamCoupled.C

\*---------------------------------------------------------------------------*/

#ifndef beamFoamCoupled_H
#define beamFoamCoupled_H

#include "sixDoFFvBeamSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace sixDoFFvBeamSolvers
{

/*---------------------------------------------------------------------------*\
                       Class beamFoamCoupled Declaration
\*---------------------------------------------------------------------------*/

class beamFoamCoupled
:
    public sixDoFFvBeamSolver
{
public:

    //- Runtime type information
    TypeName("beamFoamCoupled");


    // Constructors

        //- Construct from a dictionary and the body
        beamFoamCoupled
        (
            const dictionary& dict,
            sixDoFRigidBodyMotionFvBeam& body
        );

        //- Construct and return a clone
        virtual autoPtr<sixDoFFvBeamSolver> clone() const
        {
            return
                autoPtr<sixDoFFvBeamSolver>::NewFrom<beamFoamCoupled>
                (
                    dict_,
                    body_
                );
        }


    //- Destructor
    virtual ~beamFoamCoupled() = default;


    // Member Functions

        //- Solve the body and the beam for the fluid force and moment
        virtual void solve
        (
            bool firstIter,
            const vector& fGlobal,
            const vector& tauGlobal,
            scalar deltaT,
            scalar deltaT0
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace sixDoFFvBeamSolvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
