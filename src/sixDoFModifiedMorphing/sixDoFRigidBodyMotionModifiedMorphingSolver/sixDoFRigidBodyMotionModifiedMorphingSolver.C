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

#include "sixDoFRigidBodyMotionModifiedMorphingSolver.H"
#include "addToRunTimeSelectionTable.H"
#include "polyMesh.H"
#include "pointPatchDist.H"
#include "pointConstraints.H"
#include "boundBox.H"
#include "septernion.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sixDoFRigidBodyMotionModifiedMorphingSolver, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        sixDoFRigidBodyMotionModifiedMorphingSolver,
        dictionary
    );
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::sixDoFRigidBodyMotionModifiedMorphingSolver::cosineTransition
(
    pointScalarField& scale
)
{
    scale.primitiveFieldRef() =
        min
        (
            max
            (
                0.5
              - 0.5*cos(scale.primitiveField()*constant::mathematical::pi),
                scalar(0)
            ),
            scalar(1)
        );
}


void Foam::sixDoFRigidBodyMotionModifiedMorphingSolver::initialState
(
    point& CoR0,
    tensor& Q0
) const
{
    const sixDoFRigidBodyMotion& m = motion();

    const point c = m.transform(point::zero);
    const tensor R
    (
        tensor
        (
            m.transform(point(1, 0, 0)) - c,
            m.transform(point(0, 1, 0)) - c,
            m.transform(point(0, 0, 1)) - c
        ).T()
    );

    CoR0 = (R.T() & (m.centreOfRotation() - c));
    Q0 = (R.T() & m.orientation());
}


void Foam::sixDoFRigidBodyMotionModifiedMorphingSolver::updateXYScale()
{
    // As sixDoFRigidBodyMotionFvBeam::updateXYScale (moorFV)

    const pointField& initialPoints = points0();

    point CoR0;
    tensor Q0;
    initialState(CoR0, Q0);

    // Bounding box of the inner (rigid) region: points outside it are moved
    // to the initial centre of rotation, so they do not widen the box
    pointField points(initialPoints);

    forAll(points, pointi)
    {
        if (scale_[pointi] <= 1 - SMALL)
        {
            points[pointi] = CoR0;
        }
    }

    const boundBox box(points);
    const vector& minVal = box.min();
    const vector& maxVal = box.max();

    // Bounding box of the whole domain (suits rectangular domains)
    const boundBox domain(initialPoints);
    const vector& minDomain = domain.min();
    const vector& maxDomain = domain.max();

    // Blending lengths cut to fit inside the domain
    const scalar dx =
        min
        (
            min(minVal.x() - minDomain.x(), maxDomain.x() - maxVal.x()),
            xdist_
        );
    const scalar dy =
        min
        (
            min(minVal.y() - minDomain.y(), maxDomain.y() - maxVal.y()),
            ydist_
        );

    scalarField& xScale = xscale_.primitiveFieldRef();
    scalarField& yScale = yscale_.primitiveFieldRef();

    if (dx > SMALL)
    {
        forAll(initialPoints, pointi)
        {
            const scalar xVal = initialPoints[pointi].x();

            if (xVal >= maxVal.x())
            {
                xScale[pointi] = max(1.0 - (xVal - maxVal.x())/dx, 0.0);
            }
            else if (xVal <= minVal.x())
            {
                xScale[pointi] = max(1.0 - (minVal.x() - xVal)/dx, 0.0);
            }
            else
            {
                // Within the inner region's x range: rigid x-motion
                xScale[pointi] = 1.0;
            }
        }
    }

    if (dy > SMALL)
    {
        forAll(initialPoints, pointi)
        {
            const scalar yVal = initialPoints[pointi].y();

            if (yVal >= maxVal.y())
            {
                yScale[pointi] = max(1.0 - (yVal - maxVal.y())/dy, 0.0);
            }
            else if (yVal <= minVal.y())
            {
                yScale[pointi] = max(1.0 - (minVal.y() - yVal)/dy, 0.0);
            }
            else
            {
                yScale[pointi] = 1.0;
            }
        }

        // With an x-scale too, keep the ends in x (e.g. wave relaxation
        // zones) still
        if (dx > SMALL)
        {
            yScale *= xScale;
        }
    }

    Info<< "Modified morphing: inner region x " << minVal.x() << " to "
        << maxVal.x() << ", y " << minVal.y() << " to " << maxVal.y()
        << "; x blending length " << max(dx, scalar(0))
        << ", y blending length " << max(dy, scalar(0)) << endl;
}


Foam::tmp<Foam::pointField>
Foam::sixDoFRigidBodyMotionModifiedMorphingSolver::transform() const
{
    // As sixDoFRigidBodyMotionFvBeam::transform with x/y scales (moorFV)

    const sixDoFRigidBodyMotion& m = motion();
    const pointField& initialPoints = points0();

    point CoR0;
    tensor Q0;
    initialState(CoR0, Q0);

    const bool isXScale = xdist_ > 0;
    const bool isYScale = ydist_ > 0;

    // Translation of the body; the horizontal components handled by the
    // x/y scales are taken out of the blended (slerp) part
    const point tPoint = m.centreOfRotation() - CoR0;

    point slerpPoint(tPoint);
    if (isXScale)
    {
        slerpPoint.x() = 0;
    }
    if (isYScale)
    {
        slerpPoint.y() = 0;
    }

    const septernion s(slerpPoint, quaternion(m.orientation().T() & Q0));

    tmp<pointField> tpoints(new pointField(initialPoints));
    pointField& points = tpoints.ref();

    const scalarField& scale = scale_;
    const scalarField& xScale = xscale_;
    const scalarField& yScale = yscale_;

    forAll(points, pointi)
    {
        // Rotation and the remaining translation, blended as in the
        // standard solver
        if (scale[pointi] > SMALL)
        {
            septernion ss(s);
            if (scale[pointi] <= 1 - SMALL)
            {
                ss = slerp(septernion::I, s, scale[pointi]);
            }

            points[pointi] =
                CoR0 + ss.invTransformPoint(initialPoints[pointi] - CoR0);
        }

        // Horizontal translation
        if (isXScale && xScale[pointi] > SMALL)
        {
            points[pointi].x() += xScale[pointi]*tPoint.x();
        }
        if (isYScale && yScale[pointi] > SMALL)
        {
            points[pointi].y() += yScale[pointi]*tPoint.y();
        }
    }

    return tpoints;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionModifiedMorphingSolver::
sixDoFRigidBodyMotionModifiedMorphingSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict
)
:
    sixDoFRigidBodyMotionSolver(mesh, dict),
    xdist_(coeffDict().getOrDefault<scalar>("xDistance", -1)),
    ydist_(coeffDict().getOrDefault<scalar>("yDistance", -1)),
    scale_
    (
        IOobject
        (
            "modifiedMotionScale",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        pointMesh::New(mesh),
        dimensionedScalar(dimless, Zero)
    ),
    xscale_
    (
        IOobject
        (
            "xmotionScale",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        pointMesh::New(mesh),
        dimensionedScalar(dimless, Zero)
    ),
    yscale_
    (
        IOobject
        (
            "ymotionScale",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        pointMesh::New(mesh),
        dimensionedScalar(dimless, Zero)
    )
{
    const pointMesh& pMesh = pointMesh::New(mesh);

    // The standard solver's scale (it keeps its own copy private): 1 up to
    // innerDistance, then a cosine down to 0 at outerDistance
    const scalar di = coeffDict().get<scalar>("innerDistance");
    const scalar dOuter = coeffDict().get<scalar>("outerDistance");

    const pointPatchDist pDist
    (
        pMesh,
        mesh.boundaryMesh().patchSet(coeffDict().get<wordRes>("patches")),
        points0()
    );

    scale_.primitiveFieldRef() =
        min
        (
            max((dOuter - pDist.primitiveField())/(dOuter - di), scalar(0)),
            scalar(1)
        );

    cosineTransition(scale_);
    pointConstraints::New(pMesh).constrain(scale_);

    if (xdist_ > 0 || ydist_ > 0)
    {
        updateXYScale();

        if (xdist_ > 0)
        {
            cosineTransition(xscale_);
            pointConstraints::New(pMesh).constrain(xscale_);
            xscale_.write();
        }
        if (ydist_ > 0)
        {
            cosineTransition(yscale_);
            pointConstraints::New(pMesh).constrain(yscale_);
            yscale_.write();
        }
    }
    else
    {
        Info<< "Modified morphing: no xDistance or yDistance, so the "
            << "standard sixDoFRigidBodyMotion morphing is used" << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sixDoFRigidBodyMotionModifiedMorphingSolver::solve()
{
    // Motion update and the standard morphing
    sixDoFRigidBodyMotionSolver::solve();

    if (xdist_ > 0 || ydist_ > 0)
    {
        // Replace the standard morphing with the modified one
        pointDisplacement_.primitiveFieldRef() = transform() - points0();

        pointConstraints::New
        (
            pointDisplacement_.mesh()
        ).constrainDisplacement(pointDisplacement_);
    }
}


// ************************************************************************* //
