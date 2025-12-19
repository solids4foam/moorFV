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
#include "mathematicalConstants.H"


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

void Foam::fv::fvBeamPorosity::calculateS()
{
    // vectorField beamC;
    if (Pstream::master())
    {
        const fvMesh& beamMesh =
            mesh().time().db().parent().lookupObject<fvMesh>("beamone");

        const vectorField& beamCellCoords = beamMesh.C();   // actuator points Pa,i
        const volVectorField& refW =
            beamMesh.lookupObject<volVectorField>("refW");
        const volVectorField& W =
            beamMesh.lookupObject<volVectorField>("W");

        beamC_ = beamCellCoords + refW.internalField() + W.internalField();
    }
    Pstream::broadcast(beamC_);

    const vectorField& fluidC = mesh_.C();     // cell centres Cj

    const label nFluid = mesh_.nCells();

    sField_.setSize(nFluid);
    closestBeamCell_.setSize(nFluid);      // here: store segment index i
    closestBeamCellDist_.setSize(nFluid);  // store r

    // loop over fluid cells
    forAll(fluidC, fluidCellI)
    {
        const vector& Cj = fluidC[fluidCellI];

        scalar bestR = GREAT;
        scalar bestS = 0.0;
        label bestSeg = -1;

        // loop over beam *segments* Ei between Pa,i and Pa,i+1
        for (label i = 0; i < beamC_.size() - 1; ++i)
        {
            const vector& Pa = beamC_[i];
            const vector& Pb = beamC_[i+1];
            // segment length
            const scalar ds = mag(Pb - Pa);

            const vector Ei = Pb - Pa;
            const scalar magEi2 = magSqr(Ei);
            //skip
            if (magEi2 <= SMALL) continue;

            // parametric coordinate along the segment (can be <0 or >1)
            const scalar t = ((Cj - Pa) & Ei)/magEi2;

            // clamp to the segment, this is the normalized tangential coord s
            const scalar s = min(max(t, scalar(0)), scalar(1));

            // closest point on this segment to Cj
            const vector Pclosest = Pa + s*Ei;

            // radial distance
            const scalar r = mag(Cj - Pclosest);

            if (r < bestR)
            {
                bestR = r;
                bestS = s;
                bestSeg = i;
            }
        }

        //r
        closestBeamCellDist_[fluidCellI] = bestR;

        // segment index i
        closestBeamCell_[fluidCellI] = bestSeg;

        // s in [0,1]
        sField_[fluidCellI] = bestS;

    }
}
Foam::scalar Foam::fv::fvBeamPorosity::eta(const scalar r) const
{
    //    const scalar eps = 0.08; // epsilon ~ C * dx
    const scalar re  = r/eps_;
    return (1.0/(eps_*eps_*constant::mathematical::pi)) * exp(-re*re);
}
// Foam::tmp<Foam::volScalarField>
// Foam::fv::fvBeamPorosity::coeff(const volVectorField& U, const word& modelName) const
// {
//     auto tcoeff = tmp<volScalarField>::New
//     (
//         IOobject
//         (
//             typeName + "coeff",
//             mesh_.time().timeName(),
//             mesh_,
//             IOobject::READ_IF_PRESENT,
//             IOobject::AUTO_WRITE
//         ),
//         mesh_,
//         dimensionedScalar(dimless/dimTime, Zero)
//     );
//     auto& coeff = tcoeff.ref();
//     labelList fluidCellIDs;
//     // vectorField immersedForce;
//     label fluidRho;
//     label beamRadius;
//     if (Pstream::master())
//     {
//         fvMesh& beamMesh = const_cast<fvMesh&>(mesh().time().db().parent().lookupObject<fvMesh>("beamone"));

//         const dictionary& linkToBeamProperties = beamMesh.lookupObject<dictionary>("beamProperties");

//         fluidRho = linkToBeamProperties.subDict("coupledTotalLagNewtonRaphsonBeamCoeffs")
//         .get<dimensionedScalar>("rhoFluid").value();

//         beamRadius = linkToBeamProperties.get<dimensionedScalar>("R").value();

//         const labelList& gFluidCellIDsIO = beamMesh.lookupObject<labelIOList>("fluidCellIDs");

//         vectorIOList& immersedF = beamMesh.lookupObjectRef<vectorIOList>("immersedForce");

//         fluidCellIDs = gFluidCellIDsIO;
//         immersedForce.setSize(immersedF.size(), vector::zero);
//     }
//     else
//     {
//         fluidCellIDs.setSize(0);
//         immersedForce.setSize(0);
//         fluidRho(0);
//         beamRadius(0);
//     }
//     Pstream::broadcast(fluidCellIDs);
//     Pstream::broadcast(immersedForce);
//     Pstream::broadcast(fluidRho);
//     Pstream::broadcast(beamRadius);

//     globalIndex gI(mesh_.nCells());


//     const volScalarField& cellMarker
//     (
//         mesh_.lookupObject<volScalarField>("cellMarker")
//     );
//     const volVectorField& beamVelocity
//     (
//         mesh_.lookupObject<volVectorField>("beamVelocity")
//     );

//     if (mesh_.foundObject<volScalarField>("cellMarker"))
//     {
//         if (modelName_ == "Darcy")
//         {
//             // coeff = ((nu_/perm_) + beta_*mag(U - beamVelocity))*cellMarker;
//             coeff = (fluidRho * mag(U - beamVelocity)*beamRadius)*cellMarker;
//         }
//         // dimensions need to be fixed for fixedCoefficient model!
//         else if (modelName_ == "fixedCoefficient")
//         {
//             coeff = coeff_ * cellMarker;
//         }
//     }
//     else
//     {
//         coeff = 0.0;
//     }


//     // forAll(fluidCellIDs, beamCellI)
//     // {
//     //     const label gId = fluidCellIDs[beamCellI];
//     //     if (gId == -1)
//     //     {
//     //         continue;
//     //     }

//     //     const label owner = gI.whichProcID(gId);
//     //     // only owner should modify
//     //     if (owner != Pstream::myProcNo())
//     //     {
//     //         continue;
//     //     }

//     //     const label c = gI.toLocal(owner, gId);
//     //     if (c < 0 || c >= mesh_.nCells())
//     //     {
//     //         continue;
//     //     }
//     //     if (cellMarker[c] >= 0.001 && cellMarker[c] <= 1)
//     //     {
//     //         immersedForce[beamCellI] = coeff[c] * (U[c] - beamVelocity[c]) * mesh_.V()[c];
//     //     }
//     // }
//     // reduce(immersedForce, FieldSumOp<vector>());
//     //    Pout << "immersed Force on cell 20 = " << immersedForce[20] << endl;
//     // if (Pstream::master())
//     // {
//     //     fvMesh& beamMesh = const_cast<fvMesh&>(mesh().time().db().parent().lookupObject<fvMesh>("beamone"));

//     //     vectorIOList& immersedF =
//     //         const_cast<vectorIOList&>(beamMesh.lookupObject<vectorIOList>("immersedForce"));

//     //     if (immersedF.size() != immersedForce.size())
//     //     {
//     //         immersedF.setSize(immersedForce.size());
//     //     }
//     //     vector sumF(vector::zero);
//     //     forAll(immersedForce, i)
//     //     {
//     //         sumF += immersedForce[i];
//     //     }
//     //     Info<< "sumF on the fvOptions = " << sumF << endl;
//     //     forAll(immersedForce, i)
//     //     {
//     //         immersedF[i] = immersedForce[i];
//     //     }
//     // }

//     // vector integratedForce = vector::zero;
//     // forAll(cellMarker, cellI)
//     // {
//     //     if (cellMarker[cellI] >= 0.001 && cellMarker[cellI] <= 1)
//     //     {
//     //         integratedForce += coeff[cellI] * U[cellI] * mesh_.V()[cellI];
//     //     }
//     // }
//     // // For the integrated force, reduce before printing
//     // reduce(integratedForce, sumOp<vector>());
//     // if (forceFilePtr_.valid())
//     // {
//     //     forceFilePtr_()
//     //     << mesh_.time().timeOutputValue()
//     //     << " "
//     //     << integratedForce
//     //     << endl;
//     // }

//     coeff.correctBoundaryConditions();
//     if(mesh_.time().writeTime())
//     {
//         coeff.write();
//         // immersedF.write();
//     }
//     coeff.correctBoundaryConditions();

//     return tcoeff;
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
    fv::option(name, modelType, dict, mesh),
    forceFilePtr_(),
    eps_()
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
    Info<< "applying ALM force on fluid cells" << endl;
    // 1) Get ALM forces from beam mesh
    labelList fluidCellIDs;
    vectorField almF;

    if (Pstream::master())
    {
        const fvMesh& beamMesh =
            mesh().time().db().parent().lookupObject<fvMesh>("beamone");

        const labelList& gFluidCellIDsIO =
            beamMesh.lookupObject<labelIOList>("fluidCellIDs");

        const volVectorField& almForce =
            beamMesh.lookupObject<volVectorField>("almForce");

        almF = almForce.internalField();
        fluidCellIDs = gFluidCellIDsIO;
    }
    else
    {
        fluidCellIDs.setSize(0);
        almF.setSize(0);
    }

    Pstream::broadcast(fluidCellIDs);
    Pstream::broadcast(almF);


    // Updating s,r and segmentID's
    calculateS();
    // weights
    scalarList nodeWeights(almF.size(), 0.0);

    const scalarField& V = mesh_.V();

    // Pass 1: Accumulate weights
    forAll(closestBeamCell_, c)
    {
        const label segI = closestBeamCell_[c];

        // Check valid segment
        if (segI < 0 || segI + 1 >= almF.size())
        {
            continue;
        }

        const scalar r = closestBeamCellDist_[c];
        const scalar etaR = eta(r);

        // Optimization: skip negligible contributions
        if (etaR < SMALL)
        {
            continue;
        }

        const scalar s = sField_[c];
        const scalar cellVol = V[c];

        // Distribution logic based on Bral et al. [cite_start]Eq (9) [cite: 158]
        // Node i gets weight (1-s) * eta
        // Node i+1 gets weight (s) * eta

        nodeWeights[segI] += (1.0 - s)*etaR*cellVol;
        nodeWeights[segI + 1] += (s)*etaR*cellVol;
    }

    // Parallel Reduction: Sum weights across all processors
    // reduce(nodeWeights, sumOp<scalarList>());

    scalarField etaVals(mesh_.nCells(), 0.0);
    vectorField forceVals(mesh_.nCells(), vector::zero);

    // probably not needed now
    //    globalIndex gI(mesh_.nCells());

    // once the calculateS() has been called, apply momentum source to all
    // fluid cells, based on their s and eta
    vector ffvOption(vector::zero);
    forAll(eqn.source(), c)
    {
        const label segI = closestBeamCell_[c];
        if (segI < 0 || segI + 1 >= almF.size()) 
        {
            continue;
        }

        const scalar r = closestBeamCellDist_[c];
        scalar etaR = eta(r);
        sumEtaR += etaVals[c]*V[c];
        if (etaR < SMALL)
        {
            continue;
        }

        etaVals[c] = etaR;

        const scalar s = sField_[c];

        // We need this to convert Density (N/m) to Force (N)
        const vector& Pa = beamC_[segI];
        const vector& Pb = beamC_[segI+1];
        const scalar ds = mag(Pb - Pa);

        // Retrieve Force Density (N/m)
        const vector Fi_density   = almF[segI];
        const vector Fip1_density = almF[segI+1];

        // Normalize: (Density * Length) / Weight
        // This ensures the volume integral equals (Density * Length)
        const vector Fi_applied   = (Fi_density*ds) / (nodeWeights[segI] + VSMALL);
        const vector Fip1_applied = (Fip1_density*ds) / (nodeWeights[segI + 1] + VSMALL);

        // Interpolation and Applying the force
        const vector Sj = -etaR * ((1.0 - s) * Fi_applied + s * Fip1_applied);

        eqn.source()[c] += (Sj*V[c]);
        ffvOption += (Sj*V[c]);
    }

    Info<< "Sum of ALM forces applied to Fluid (fvOptions) = " << ffvOption << endl;
    if (mesh_.time().writeTime())
    {
        volScalarField etaField
        (
            IOobject
            (
                "eta",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimless/dimArea, 0.0),
            "zeroGradient"
        );

        //        etaField.internalField() = etaVals;
        forAll(etaVals, c)
        {
            etaField.internalFieldRef()[c] = etaVals[c];
        }
        etaField.correctBoundaryConditions();
        etaField.write();
        volVectorField forceField
        (
            IOobject
            (
                "ffield",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedVector("zero", dimForce, vector::zero),
            "zeroGradient"
        );

        //        etaField.internalField() = etaVals;
        forAll(forceField, c)
        {
            forceField.internalFieldRef()[c] = forceVals[c];
        }
        forceField.correctBoundaryConditions();
        forceField.write();
}
}


void Foam::fv::fvBeamPorosity::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    const volVectorField& U = eqn.psi();

    // fvMatrix<vector> mangrovesEqn
    // (
    //   - fvm::Sp(rho*coeff(U, modelName_), U)
    // );

    // Contributions are added to RHS of momentum equation
    // eqn += mangrovesEqn;
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
        coeffs_.readEntry("epsilon", eps_);
        return true;
    }
    return false;
}


// ************************************************************************* //
