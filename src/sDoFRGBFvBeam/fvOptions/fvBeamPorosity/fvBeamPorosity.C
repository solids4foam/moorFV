/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 IH-Cantabria
    Copyright (C) 2017-2021 OpenCFD Ltd.
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

#include "fvBeamPorosity.H"
#include "mathematicalConstants.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "fvmDdt.H"
#include "fvmSup.H"
#include "addToRunTimeSelectionTable.H"
#include "samplingFluid.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(fvBeamPorosity, 0);
    addToRunTimeSelectionTable
    (
        option,
        fvBeamPorosity,
        dictionary
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::fv::fvBeamPorosity::dragCoeff(const volVectorField& U) const
{
    auto tdragCoeff = tmp<volScalarField>::New
    (
        IOobject
        (
            typeName + ":dragCoeff",
            mesh_.time().timeName(),
            mesh_.time(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/dimTime, Zero)
    );
    auto& dragCoeff = tdragCoeff.ref();
    const dictionary& linkToBeamProperties = mesh_.time().db().parent()
        .lookupObject<dictionary>("beamProperties");
    const dimensionedScalar& fluidRho =
        linkToBeamProperties.subDict("coupledTotalLagNewtonRaphsonBeamCoeffs")
            .get<dimensionedScalar>("rhoFluid");
    const dimensionedScalar& beamRadius = linkToBeamProperties
        .get<dimensionedScalar>("R");
    const dimensionedScalar& beamLength = linkToBeamProperties
        .get<dimensionedScalar>("L");
    const volScalarField& cellMarker
    (
        mesh_.lookupObject<volScalarField>("cellMarker")
    );
    forAll(mesh_.C(),celli)
    {
        if (mesh_.foundObject<volScalarField>("cellMarker"))
        {
            dragCoeff[celli] = 0.5 * cellMarker[celli] * fluidRho.value()
              * (2 * beamRadius.value()) * beamLength.value() * mag(U[celli]);
        }
        else
        {
            dragCoeff[celli] = 0;
        }
    }
    dragCoeff.correctBoundaryConditions();
    return tdragCoeff;
}


Foam::tmp<Foam::volScalarField>
Foam::fv::fvBeamPorosity::inertiaCoeff() const
{
    auto tinertiaCoeff = tmp<volScalarField>::New
    (
        IOobject
        (
            typeName + ":inertiaCoeff",
            mesh_.time().timeName(),
            mesh_.time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    );
    auto& inertiaCoeff = tinertiaCoeff.ref();

    const scalar pi = constant::mathematical::pi;
    const dictionary& linkToBeamProperties = mesh_.time().db().parent()
        .lookupObject<dictionary>("beamProperties");
    const dimensionedScalar& fluidRho =
        linkToBeamProperties.subDict("coupledTotalLagNewtonRaphsonBeamCoeffs")
            .get<dimensionedScalar>("rhoFluid");
    const dimensionedScalar& beamRadius = linkToBeamProperties
        .get<dimensionedScalar>("R");
    const dimensionedScalar& beamLength = linkToBeamProperties
        .get<dimensionedScalar>("L");
    const scalar& cm = linkToBeamProperties.getOrDefault<scalar>("CMn", 1.0);
    const volScalarField& cellMarker
    (
        mesh_.lookupObject<volScalarField>("cellMarker")
    );

    forAll(mesh_.C(),celli)
    {
        if (mesh_.foundObject<volScalarField>("cellMarker"))
        {
            inertiaCoeff[celli] = (cm+1) * fluidRho.value() * cellMarker[celli]
              * beamLength.value() * pow((2 * beamRadius.value()),2) * pi;
        }
        else
        {
            inertiaCoeff[celli] = 0;
        }
    }
    inertiaCoeff.correctBoundaryConditions();

    return tinertiaCoeff;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::fvBeamPorosity::fvBeamPorosity
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fv::option(name, modelType, dict, mesh)
    // active_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::fvBeamPorosity::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    const volVectorField& U = eqn.psi();

    fvMatrix<vector> mangrovesEqn
    (
      - fvm::Sp(dragCoeff(U), U)
      - inertiaCoeff()*fvm::ddt(U)
    );
    // Contributions are added to RHS of momentum equation
    eqn += mangrovesEqn;
}


void Foam::fv::fvBeamPorosity::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    const volVectorField& U = eqn.psi();

    fvMatrix<vector> mangrovesEqn
    (
      - fvm::Sp(rho*dragCoeff(U), U)
      - rho*inertiaCoeff()*fvm::ddt(U)
    );

    // Contributions are added to RHS of momentum equation
    eqn += mangrovesEqn;
}


bool Foam::fv::fvBeamPorosity::read(const dictionary& dict)
{
    if (fv::option::read(dict))
    {
            if (!coeffs_.readIfPresent("UNames", fieldNames_))
            {
                fieldNames_.resize(1);
                fieldNames_.first() = coeffs_.getOrDefault<word>("U", "U");
            }
        fv::option::resetApplied();
        return true;
    }
    return false;
}


// ************************************************************************* //
