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

\*---------------------------------------------------------------------------*/

#include "sixDoFRigidBodyMotionFvBeamSolver.H"
#include "addToRunTimeSelectionTable.H"
#include "polyMesh.H"
#include "pointPatchDist.H"
#include "pointConstraints.H"
#include "uniformDimensionedFields.H"
#include "forces.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sixDoFRigidBodyMotionFvBeamSolver, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        sixDoFRigidBodyMotionFvBeamSolver,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionFvBeamSolver::sixDoFRigidBodyMotionFvBeamSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict
)
:
    displacementMotionSolver(mesh, dict, typeName),
    motion_
    (
        coeffDict(),
        IOobject
        (
            "sixDoFRigidBodyMotionFvBeamState",
            mesh.time().timeName(),
            "uniform",
            mesh
        ).typeHeaderOk<IOdictionary>(true)
      ? IOdictionary
        (
            IOobject
            (
                "sixDoFRigidBodyMotionFvBeamState",
                mesh.time().timeName(),
                "uniform",
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE,
                false
            )
        )
      : coeffDict(),
        mesh.time()
    ),
    patches_(coeffDict().get<wordRes>("patches")),
    patchSet_(mesh.boundaryMesh().patchSet(patches_)),
    di_(coeffDict().get<scalar>("innerDistance")),
    do_(coeffDict().get<scalar>("outerDistance")),
    test_(coeffDict().getOrDefault("test", false)),
    rhoInf_(1.0),
    rhoName_(coeffDict().getOrDefault<word>("rho", "rho")),
    scale_
    (
        IOobject
        (
            "motionScale",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        pointMesh::New(mesh),
        dimensionedScalar(dimless, Zero)
    ),
    curTimeIndex_(-1),
    CofGvelocity_(coeffDict().getOrDefault<word>("CofGvelocity", "none"))
{
    if (rhoName_ == "rhoInf")
    {
        coeffDict().readEntry("rhoInf", rhoInf_);
    }

    // Calculate scaling factor everywhere
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    {
        const pointMesh& pMesh = pointMesh::New(mesh);

        pointPatchDist pDist(pMesh, patchSet_, points0());

        if (mag(do_ - di_) < SMALL || do_ < di_)
        {
            FatalErrorInFunction 
                << "Invalid values: Outer diameter (do_) and inner diameter (di_) are too close or incorrectly defined.\n"
                << "Ensure that do_ > di_ with sufficient difference (greater than SMALL).\n"
                << "Current values: do_ = " << do_ << ", di_ = " << di_ 
                << abort(FatalError);
        }

        // Scaling: 1 up to di then linear down to 0 at do away from patches
        scale_.primitiveFieldRef() =
            min
            (
                max
                (
                    (do_ - pDist.primitiveField())/(do_ - di_),
                    scalar(0)
                ),
                scalar(1)
            );

        // Convert the scale function to a cosine
        scale_.primitiveFieldRef() =
            min
            (
                max
                (
                    0.5
                  - 0.5
                   *cos(scale_.primitiveField()
                   *Foam::constant::mathematical::pi),
                    scalar(0)
                ),
                scalar(1)
            );

        pointConstraints::New(pMesh).constrain(scale_);
        scale_.write();
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::sixDoFRigidBodyMotionFvBeam&
Foam::sixDoFRigidBodyMotionFvBeamSolver::motion() const
{
    return motion_;
}


Foam::tmp<Foam::pointField>
Foam::sixDoFRigidBodyMotionFvBeamSolver::curPoints() const
{
    tmp<pointField> newPoints
    (
        points0() + pointDisplacement_.primitiveField()
    );

    if (!moveAllCells())
    {
        tmp<pointField> ttransformedPts(new pointField(mesh().points()));
        pointField& transformedPts = ttransformedPts.ref();

        UIndirectList<point>(transformedPts, pointIDs()) =
            pointField(newPoints.ref(), pointIDs());

        return ttransformedPts;
    }

    return newPoints;
}


void Foam::sixDoFRigidBodyMotionFvBeamSolver::solve()
{
    const Time& t = mesh().time();

    if (mesh().nPoints() != points0().size())
    {
        FatalErrorInFunction
            << "The number of points in the mesh seems to have changed." << endl
            << "In constant/polyMesh there are " << points0().size()
            << " points; in the current mesh there are " << mesh().nPoints()
            << " points." << exit(FatalError);
    }

    // Store the motion state at the beginning of the time-step
    bool firstIter = false;
    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        motion_.newTime();
        curTimeIndex_ = this->db().time().timeIndex();
        firstIter = true;
    }

    dimensionedVector g("g", dimAcceleration, Zero);

    if (db().time().foundObject<uniformDimensionedVectorField>("g"))
    {
        g = db().time().lookupObject<uniformDimensionedVectorField>("g");
    }
    else
    {
        coeffDict().readIfPresent("g", g);
    }

    // const scalar ramp = min(max((this->db().time().value() - 5)/10, 0), 1);
    const scalar ramp = 1.0;

    if (test_)
    {
        motion_.update
        (
            firstIter,
            ramp*(motion_.mass()*g.value()),
            ramp*(motion_.mass()*(motion_.momentArm() ^ g.value())),
            t.deltaTValue(),
            t.deltaT0Value()
        );
    }
    else
    {
        dictionary forcesDict;

        forcesDict.add("type", functionObjects::forces::typeName);
        forcesDict.add("patches", patches_);
        forcesDict.add("rhoInf", rhoInf_);
        forcesDict.add("rho", rhoName_);
        forcesDict.add("CofR", motion_.centreOfRotation());

        functionObjects::forces f("forces", db(), forcesDict);

        f.calcForcesMoments();

        motion_.update
        (
            firstIter,
            ramp*(f.forceEff() + motion_.mass()*g.value()),
            ramp
           *(
               f.momentEff()
             + motion_.mass()*(motion_.momentArm() ^ g.value())
            ),
            t.deltaTValue(),
            t.deltaT0Value()
        );

        if (CofGvelocity_ != "none")
        {
            if
            (
                db().time().foundObject<uniformDimensionedVectorField>
                (
                    CofGvelocity_
                )
            )
            {
                uniformDimensionedVectorField& vel =
                    db().time().lookupObjectRef<uniformDimensionedVectorField>
                    (
                        CofGvelocity_
                    );
                vel = motion_.v();
            }
        }
    }

    // Update the displacements
    pointDisplacement_.primitiveFieldRef() =
        motion_.transform(points0(), scale_) - points0();

    // Displacement has changed. Update boundary conditions
    pointConstraints::New
    (
        pointDisplacement_.mesh()
    ).constrainDisplacement(pointDisplacement_);
}


bool Foam::sixDoFRigidBodyMotionFvBeamSolver::writeObject
(
    IOstreamOption streamOpt,
    const bool valid
) const
{
    IOdictionary dict
    (
        IOobject
        (
            "sixDoFRigidBodyMotionFvBeamState",
            mesh().time().timeName(),
            "uniform",
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    motion_.state().write(dict);
    return dict.regIOobject::write();
}


bool Foam::sixDoFRigidBodyMotionFvBeamSolver::read()
{
    if (displacementMotionSolver::read())
    {
        motion_.read(coeffDict(), db().time());

        return true;
    }

    return false;
}


// ************************************************************************* //
