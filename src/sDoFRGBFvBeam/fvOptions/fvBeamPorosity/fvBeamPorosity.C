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
#include "couplingHelperFunctions.H"


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
    label celli;
    auto& dragCoeff = tdragCoeff.ref();
    
    Info << "available objects = " << U.db().names() << endl;
    if (mesh_.foundObject<volScalarField>("fluidCellMarker"))
    {
        Info << " im here " << endl;
        const volScalarField& gg(mesh_.lookupObject<volScalarField>("fluidCellMarker"));
        forAll(mesh_.C(),celli)
        {
            dragCoeff[celli] = 2000*gg[celli];
        }
        if (mesh_.time().writeTime())
        {
            dragCoeff.write();
        }
    }
    else
    {
        dragCoeff[celli] = 0;
    }
    

    return tdragCoeff;
}


// Foam::tmp<Foam::volScalarField>
// Foam::fv::fvBeamPorosity::inertiaCoeff() const
// {
//     auto tinertiaCoeff = tmp<volScalarField>::New
//     (
//         IOobject
//         (
//             typeName + ":inertiaCoeff",
//             mesh_.time().timeName(),
//             mesh_.time(),
//             IOobject::NO_READ,
//             IOobject::NO_WRITE
//         ),
//         mesh_,
//         dimensionedScalar(dimless, Zero)
//     );
//     auto& inertiaCoeff = tinertiaCoeff.ref();

//     const scalar pi = constant::mathematical::pi;

//     forAll(zoneIDs_, i)
//     {
//         const scalar a = aZone_[i];
//         const scalar N = NZone_[i];
//         const scalar Cm = CmZone_[i];

//         const labelList& zones = zoneIDs_[i];

//         for (label zonei : zones)
//         {
//             const cellZone& cz = mesh_.cellZones()[zonei];

//             for (label celli : cz)
//             {
//                 inertiaCoeff[celli] = 0.25*(Cm+1)*pi*a*a*N;
//             }
//         }
//     }

//     inertiaCoeff.correctBoundaryConditions();

//     return tinertiaCoeff;
// }


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
      //- inertiaCoeff()*fvm::ddt(U)
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
      //- rho*inertiaCoeff()*fvm::ddt(U)
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

    //     // Create the Mangroves models - 1 per region
    //     const dictionary& regionsDict(coeffs_.subDict("regions"));
    //     const wordList regionNames(regionsDict.toc());
    //     aZone_.setSize(regionNames.size(), 1);
    //     NZone_.setSize(regionNames.size(), 1);
    //     CmZone_.setSize(regionNames.size(), 1);
    //     CdZone_.setSize(regionNames.size(), 1);
    //     zoneIDs_.setSize(regionNames.size());

    //     forAll(zoneIDs_, i)
    //     {
    //         const word& regionName = regionNames[i];
    //         const dictionary& modelDict = regionsDict.subDict(regionName);

    //         const word zoneName(modelDict.get<word>("cellZone"));

    //         zoneIDs_[i] = mesh_.cellZones().indices(zoneName);
    //         if (zoneIDs_[i].empty())
    //         {
    //             FatalErrorInFunction
    //                 << "Unable to find cellZone " << zoneName << nl
    //                 << "Valid cellZones are:" << mesh_.cellZones().names()
    //                 << exit(FatalError);
    //         }

    //         modelDict.readEntry("a", aZone_[i]);
    //         modelDict.readEntry("N", NZone_[i]);
    //         modelDict.readEntry("Cm", CmZone_[i]);
    //         modelDict.readEntry("Cd", CdZone_[i]);
    //     }

        return true;
    }

    return false;
}


// ************************************************************************* //
