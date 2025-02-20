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
Foam::fv::fvBeamPorosity::coeff(const volVectorField& U, const word& modelName) const
{
    auto tcoeff = tmp<volScalarField>::New
    (
        IOobject
        (
            typeName + "coeff",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/dimTime, Zero)
    );
    auto& coeff = tcoeff.ref();

    const volScalarField& cellMarker
    (
        mesh_.lookupObject<volScalarField>("cellMarker")
    );

    forAll(mesh_.C(),celli)
    {
        if (mesh_.foundObject<volScalarField>("cellMarker"))
        {
            if (modelName_ == "DarcyLike")
            {
                coeff[celli] = (nu_ / perm_) * cellMarker[celli];
            }
            if (modelName_ == "SmagorinskyLike")
            {
                coeff[celli] = rho_ * pFactor_ * cellMarker[celli] * pow(mag(U[celli]),exponent_);
            }
        }
        else
        {
            coeff[celli] = 0;
        }
    }
    vector integratedForce = vector::zero;
    forAll(cellMarker, cellI)
    {
        if (cellMarker[cellI] >= 0.001 && cellMarker[cellI] <= 1)
        {
            integratedForce += coeff[cellI] * U[cellI] * mesh_.V()[cellI];
        }
    }
    if (forceFilePtr_.valid())
    {
        forceFilePtr_()
        << mesh_.time().timeOutputValue()
        << " "
        << integratedForce
        << endl;
    }

    coeff.correctBoundaryConditions();
    if(mesh_.time().writeTime())
    {
        coeff.write();
    }
    coeff.correctBoundaryConditions();

    return tcoeff;
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
    fv::option(name, modelType, dict, mesh),
    forceFilePtr_(),
    modelName_(coeffs_.get<word>("modelName")),
    perm_(),
    nu_(),
    rho_(),
    pFactor_(),
    exponent_()
    // active_()
{
    read(dict);
    fileName postProcDir;
    word startTimeName =
            mesh.time().timeName(mesh.time().startTime().value());
    postProcDir = mesh.time().path()/"postProcessing"/startTimeName;
    Info <<" postProcDir =  " << postProcDir << endl;


    mkDir(postProcDir);
    forceFilePtr_.reset
    (
        new OFstream
        (
            postProcDir/"dragCoeffIntegration.dat"
        )
    );
    if (forceFilePtr_.valid())
    {
        forceFilePtr_()
            << "# Time"
            << " "
            << "force"
            << endl;
    }
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
      - fvm::Sp(coeff(U, modelName_), U)
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
      - fvm::Sp(rho*coeff(U, modelName_), U)
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
        Info << "modelName = " << modelName_ << endl;
        if (modelName_ == "DarcyLike")
        {
            coeffs_.readEntry("k", perm_);
            coeffs_.readEntry("nu", nu_);
        }
        else if (modelName_ == "SmagorinskyLike")
        {
            coeffs_.readEntry("rho", rho_);
            coeffs_.readEntry("penalizationFactor", pFactor_);
            coeffs_.readEntry("exponent", exponent_);
        }
        else
        {
            FatalError
            (
                "Foam::fv::fvBeamPorosity::read"
            )
            << "Unknown model type " << modelName_ << nl
            << "Valid model types are: DarcyLike, smagorinskyLike" << nl
            << abort(FatalError);
        }

        return true;
    }
    return false;
}


// ************************************************************************* //
