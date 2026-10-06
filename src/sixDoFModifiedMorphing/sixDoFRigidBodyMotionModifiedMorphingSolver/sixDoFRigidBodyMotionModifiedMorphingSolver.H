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
    Foam::sixDoFRigidBodyMotionModifiedMorphingSolver

Description
    The standard sixDoFRigidBodyMotion mesh-motion solver with moorFV's
    modified mesh morphing (as in sixDoFRigidBodyMotionFvBeamSolver), so that
    bodies moved by the standard sixDoFRigidBodyMotion (and its restraints,
    e.g. the foamMooring moorDynR2 restraint) can use it.

    The rotation and the vertical translation are blended over
    innerDistance/outerDistance as in the standard solver. The horizontal
    translation is instead applied rigidly to every point within the
    x (and y) range of the inner region, and blended to zero over xDistance
    (yDistance) beyond it, with a cosine transition, through the whole depth
    of the domain. The blending lengths are cut to fit inside the domain.
    This spreads large surge (sway) motions over many more cells than the
    standard morphing.

    The motion, forces, restraints, constraints, state file and report are
    those of the standard solver, which reads sixDoFRigidBodyMotionCoeffs.

Usage
    \verbatim
    motionSolverLibs    (sixDoFRigidBodyMotion sixDoFModifiedMorphing);

    motionSolver        sixDoFRigidBodyMotionModifiedMorphing;

    sixDoFRigidBodyMotionCoeffs
    {
        // As for sixDoFRigidBodyMotion, plus:
        xDistance   0.3;    // x blending length; < 0 or absent: not used
        yDistance   0.15;   // y blending length; < 0 or absent: not used
    }
    \endverbatim

    Without xDistance and yDistance the motion is that of the standard
    solver.

SourceFiles
    sixDoFRigidBodyMotionModifiedMorphingSolver.C

\*---------------------------------------------------------------------------*/

#ifndef sixDoFRigidBodyMotionModifiedMorphingSolver_H
#define sixDoFRigidBodyMotionModifiedMorphingSolver_H

#include "sixDoFRigidBodyMotionSolver.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class sixDoFRigidBodyMotionModifiedMorphingSolver Declaration
\*---------------------------------------------------------------------------*/

class sixDoFRigidBodyMotionModifiedMorphingSolver
:
    public sixDoFRigidBodyMotionSolver
{
    // Private Data

        //- Horizontal blending lengths (< 0: not used)
        const scalar xdist_;
        const scalar ydist_;

        //- Rotation and vertical-translation scale: 1 within innerDistance,
        //  0 beyond outerDistance (as the standard solver's motionScale)
        pointScalarField scale_;

        //- Horizontal-translation scales
        pointScalarField xscale_;
        pointScalarField yscale_;


    // Private Member Functions

        //- Turn a linear 0-1 transition into a cosine one
        static void cosineTransition(pointScalarField& scale);

        //- Initial centre of rotation and orientation of the motion (not
        //  public in sixDoFRigidBodyMotion), recovered from its public
        //  rigid-body transform: transform(p) = c + R p with R = Q Q0^T
        void initialState(point& CoR0, tensor& Q0) const;

        //- Set the horizontal-translation scales from the initial points
        void updateXYScale();

        //- Points moved by the current motion state with the modified
        //  morphing
        tmp<pointField> transform() const;

        //- No copy construct
        sixDoFRigidBodyMotionModifiedMorphingSolver
        (
            const sixDoFRigidBodyMotionModifiedMorphingSolver&
        ) = delete;

        //- No copy assignment
        void operator=
        (
            const sixDoFRigidBodyMotionModifiedMorphingSolver&
        ) = delete;


public:

    //- Runtime type information
    TypeName("sixDoFRigidBodyMotionModifiedMorphing");


    // Constructors

        //- Construct from polyMesh and IOdictionary
        sixDoFRigidBodyMotionModifiedMorphingSolver
        (
            const polyMesh&,
            const IOdictionary& dict
        );


    //- Destructor
    virtual ~sixDoFRigidBodyMotionModifiedMorphingSolver() = default;


    // Member Functions

        //- Solve for the motion, then morph the mesh with the modified
        //  morphing
        virtual void solve();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
