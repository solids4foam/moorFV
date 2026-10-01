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

\*---------------------------------------------------------------------------*/

#include "beamFoamCoupled.H"
#include "finiteVolumeBeam.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sixDoFFvBeamSolvers
{
    defineTypeNameAndDebug(beamFoamCoupled, 0);
    addToRunTimeSelectionTable(sixDoFFvBeamSolver, beamFoamCoupled, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDoFFvBeamSolvers::beamFoamCoupled::beamFoamCoupled
(
    const dictionary& dict,
    sixDoFRigidBodyMotionFvBeam& body
)
:
    sixDoFFvBeamSolver(dict, body)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sixDoFFvBeamSolvers::beamFoamCoupled::solve
(
    bool firstIter,
    const vector& fGlobal,
    const vector& tauGlobal,
    scalar deltaT,
    scalar deltaT0
)
{
    typedef sixDoFRigidBodyMotionFvBeamRestraints::finiteVolumeBeam
        finiteVolumeBeam;

    if (constraints().size())
    {
        FatalErrorInFunction
            << "beamFoamCoupled does not support motion constraints"
            << abort(FatalError);
    }

    // Applied load: fluid force and moment with the weight, plus any other
    // restraints at the current state (explicit). Moments about the centre
    // of rotation, which is the centre of mass
    const point CoR = centreOfRotation();

    vector force = fGlobal;
    vector moment = tauGlobal;

    const finiteVolumeBeam* beamRestraint = nullptr;

    forAll(restraints(), rI)
    {
        const sixDoFRigidBodyMotionFvBeamRestraint& restraint =
            restraints()[rI];

        if (isA<finiteVolumeBeam>(restraint))
        {
            if (beamRestraint)
            {
                FatalErrorInFunction
                    << "beamFoamCoupled supports one finiteVolumeBeam "
                    << "restraint" << abort(FatalError);
            }

            beamRestraint = &refCast<const finiteVolumeBeam>(restraint);
            continue;
        }

        point rP = Zero;
        vector rF = Zero;
        vector rM = Zero;

        restraint.restrain(body_, rP, rF, rM);

        force += rF;
        moment += rM + ((rP - CoR) ^ rF);
    }

    if (!beamRestraint)
    {
        FatalErrorInFunction
            << "beamFoamCoupled needs a finiteVolumeBeam restraint"
            << abort(FatalError);
    }

    const RigidBodyEndState state =
        beamRestraint->solveMonolithic(body_, force, moment, dict_);

    // Accept beamFoam's state as it is. The motion stores the angular
    // momentum and torque in body axes
    const tensor J
    (
        body_.momentOfInertia().xx(), 0, 0,
        0, body_.momentOfInertia().yy(), 0,
        0, 0, body_.momentOfInertia().zz()
    );

    centreOfRotation() = initialCentreOfRotation() + state.displacement;
    Q() = state.orientation;
    v() = state.velocity;
    a() = state.acceleration;
    pi() = (J & (state.orientation.T() & state.angularVelocity));
    tau() = (state.orientation.T() & (moment + state.beamTorque));

    Info<< "beamFoamCoupled: centre of rotation " << centreOfRotation()
        << ", velocity " << v()
        << ", angular velocity " << state.angularVelocity
        << ", beam force " << state.beamForce
        << ", Newton iterations " << state.iterations << endl;
}


// ************************************************************************* //
