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
#include "vectorIOList.H"
#include "FieldSumOp.H"
#include "fvcSup.H"

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
    labelList fluidCellIDs;
    vectorField immersedForce;
    if (Pstream::master())
    {
        fvMesh& beamMesh = const_cast<fvMesh&>(mesh().time().db().parent().lookupObject<fvMesh>("beamone"));

        const labelList& gFluidCellIDsIO = beamMesh.lookupObject<labelIOList>("fluidCellIDs");

        vectorIOList& immersedF = beamMesh.lookupObjectRef<vectorIOList>("immersedForce");

        fluidCellIDs = gFluidCellIDsIO;
        immersedForce.setSize(immersedF.size(), vector::zero);
    }
    else
    {
        fluidCellIDs.setSize(0);
        immersedForce.setSize(0);
    }
    Pstream::broadcast(fluidCellIDs);
    Pstream::broadcast(immersedForce);

    globalIndex gI(mesh_.nCells());


    const volScalarField& cellMarker
    (
        mesh_.lookupObject<volScalarField>("cellMarker")
    );
    const volVectorField& beamVelocity
    (
        mesh_.lookupObject<volVectorField>("beamVelocity")
    );

    if (mesh_.foundObject<volScalarField>("cellMarker"))
    {
        if (modelName_ == "Darcy")
        {
            coeff = (nu_/perm_) + beta_*mag(U - beamVelocity)  * cellMarker;
        }
        // dimensions need to be fixed for fixedCoefficient model!
        else if (modelName_ == "fixedCoefficient")
        {
            coeff = coeff_ * cellMarker;
        }
    }
    else
    {
        coeff = 0.0;
    }


    forAll(fluidCellIDs, beamCellI)
    {
        const label gId = fluidCellIDs[beamCellI];
        if (gId == -1)
        {
            continue;
        }

        const label owner = gI.whichProcID(gId);
        // only owner should modify
        if (owner != Pstream::myProcNo())
        {
            continue;
        }

        const label c = gI.toLocal(owner, gId);
        if (c < 0 || c >= mesh_.nCells())
        {
            continue;
        }
        if (cellMarker[c] >= 0.001 && cellMarker[c] <= 1)
        {
            immersedForce[beamCellI] = coeff[c] * (U[c] - beamVelocity[c]) * mesh_.V()[c];
        }
    }
    reduce(immersedForce, FieldSumOp<vector>());
    //    Pout << "immersed Force on cell 20 = " << immersedForce[20] << endl;
    if (Pstream::master())
    {
        fvMesh& beamMesh = const_cast<fvMesh&>(mesh().time().db().parent().lookupObject<fvMesh>("beamone"));

        vectorIOList& immersedF =
            const_cast<vectorIOList&>(beamMesh.lookupObject<vectorIOList>("immersedForce"));

        if (immersedF.size() != immersedForce.size())
        {
            immersedF.setSize(immersedForce.size());
        }
        vector sumF(vector::zero);
        forAll(immersedForce, i)
        {
            sumF += immersedForce[i];
        }
        Info<< "sumF on the fvOptions = " << sumF << endl;
        forAll(immersedForce, i)
        {
            immersedF[i] = immersedForce[i];
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
    // For the integrated force, reduce before printing
    reduce(integratedForce, sumOp<vector>());
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
        // immersedF.write();
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
    beta_(),
    coeff_()
    // active_()
{
    read(dict);

    if (Pstream::master())
    {
        fileName postProcDir;

        word startTimeName =
                mesh.time().timeName(mesh.time().startTime().value());

        if (Pstream::parRun())
        {
            // Put in undecomposed case (Note: gives problems for
            // distributed data running)
            postProcDir = mesh.time().path()/".."/"postProcessing"/startTimeName;
        }
        else
        {
            postProcDir = mesh.time().path()/"postProcessing"/startTimeName;
        }

        mkDir(postProcDir);
        forceFilePtr_.reset
        (
            new OFstream
            (
                postProcDir/"integratedForce.dat"
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
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::fvBeamPorosity::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
    {
    const volVectorField& U = eqn.psi();
    const volVectorField& Ubeam = U.mesh().lookupObject<volVectorField>("beamVelocity");


    fvMatrix<vector> mangrovesEqn
    (
        - fvm::Sp(coeff(U, modelName_), U)
    );
    // Contributions are added to RHS of momentum equation
    eqn += mangrovesEqn;
    eqn += Ubeam*coeff(U, modelName_);
    // eqn += fvc::SuSp(coeff(U, modelName_), Ubeam);
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
        if (modelName_ == "Darcy")
        {
            coeffs_.readEntry("k", perm_);
            coeffs_.readEntry("nu", nu_);
            coeffs_.readEntry("beta", beta_);
        }
        else if (modelName_ == "fixedCoefficient")
        {
            coeffs_.readEntry("coefficient", coeff_);
        }
        else
        {
            FatalError
            (
                "Foam::fv::fvBeamPorosity::read"
            )
            << "Unknown model type " << modelName_ << nl
            << "Valid model types are: Darcy,  fixedCoefficient" << nl
            << abort(FatalError);
        }

        return true;
    }
    return false;
}


// ************************************************************************* //
