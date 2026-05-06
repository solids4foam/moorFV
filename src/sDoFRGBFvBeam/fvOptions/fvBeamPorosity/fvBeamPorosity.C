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

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //
//------------------------------------------------------------------------------
/*
 * Compute geometric mapping between fluid cells and beam segments.
 *
 * For each fluid cell centre:
 *  - Finds the closest beam segment
 *  - Computes the shortest radial distance to the segment
 *  - Computes the local interpolation coordinate (s in [0,1])
 *
 * Outputs:
 *  - closestBeamCell_  : index of nearest beam segment
 *  - closestBeamCellDist_ : radial distance from beam centreline
 *  - sField_           : interpolation parameter along segment
 *
 * Complexity: O(N_fluid * N_beamSegments)
 */
//------------------------------------------------------------------------------

void fvBeamPorosity::calculateS()
{
    if (Pstream::master())
    {
        const fvMesh& beamMesh =
            mesh().time().db().parent().lookupObject<fvMesh>(beamName_);

        // Actuator points Pa,i
        const vectorField& beamCellCoords = beamMesh.C();
        const volVectorField& refW =
            beamMesh.lookupObject<volVectorField>("refW");
        const volVectorField& W =
            beamMesh.lookupObject<volVectorField>("W");

        beamC_ = beamCellCoords + refW.internalField() + W.internalField();
    }

    Pstream::broadcast(beamC_);
    // cell centres Cj
    const vectorField& fluidC = mesh_.C();

    const label nFluid = mesh_.nCells();

    sField_.setSize(nFluid);
    // here: store segment index i
    closestBeamCell_.setSize(nFluid);
    // store r
    closestBeamCellDist_.setSize(nFluid);

    // boundBox beamBb(beamC_, false);
    // beamBb.inflate(3.0*eps_);

    // loop over fluid cells
    forAll(fluidC, fluidCellI)
    {
        const vector& Cj = fluidC[fluidCellI];

        // if (!beamBb.contains(Cj))
        // {
        //     closestBeamCellDist_[fluidCellI] = GREAT;
        //     closestBeamCell_[fluidCellI] = -1;
        //     sField_[fluidCellI] = 0.0;
        //     continue;
        // }

        scalar bestR = GREAT;
        scalar bestS = 0.0;
        label bestSeg = -1;

        // loop over beam *segments* Ei between Pa,i and Pa,i+1
        for (label i = 0; i < beamC_.size() - 1; ++i)
        {
            const vector& Pa = beamC_[i];
            const vector& Pb = beamC_[i + 1];
            const vector Ei = Pb - Pa;
            const scalar magEi2 = magSqr(Ei);
            //skip
            if (magEi2 <= SMALL)
            {
                continue;
            }

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

// Gaussian kernel function
scalar fvBeamPorosity::eta(const scalar r) const
{
    const scalar re = r/eps_;
    return (1.0/(eps_*eps_*constant::mathematical::pi))*exp(-re*re);
}

// Apply beam force to fluid cell
void fvBeamPorosity::applyBeamForce(fvMatrix<vector>& eqn)
{
    //    const volVectorField& U = eqn.psi();
    Info<< "applying ALM force on fluid cells" << endl;

    // 1) Get ALM forces from beam mesh
    labelList fluidCellIDs;
    vectorField almF;

    if (Pstream::master())
    {
        const fvMesh& beamMesh =
            mesh().time().db().parent().lookupObject<fvMesh>(beamName_);

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

        // if (r > 3.0*eps_)
        // {
        //     continue;
        // }

        const scalar etaR = eta(r);

        // Optimization: skip negligible contributions
        if (etaR < SMALL)
        {
            continue;
        }

        const scalar s = sField_[c];
        const scalar cellVol = V[c];

        // Distribution logic based on Bral et al.
        // Node i gets weight (1-s)*eta
        // Node i+1 gets weight (s)*eta

        nodeWeights[segI] += (1.0 - s)*etaR*cellVol;
        nodeWeights[segI + 1] += s*etaR*cellVol;
    }

    // Parallel Reduction: Sum weights across all processors
    reduce(nodeWeights, sumOp<scalarList>());

    scalarField etaVals(mesh_.nCells(), 0.0);
    vectorField forceVals(mesh_.nCells(), vector::zero);

    // probably not needed now
    // globalIndex gI(mesh_.nCells());

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

        // if (r > 3.0*eps_)
        // {
        //     continue;
        // }

        scalar etaR = eta(r);
        if (etaR < SMALL)
        {
            continue;
        }

        etaVals[c] = etaR;

        const scalar s = sField_[c];

        // We need this to convert Density (N/m) to Force (N)
        const vector& Pa = beamC_[segI];
        const vector& Pb = beamC_[segI + 1];
        const scalar ds = mag(Pb - Pa);

        // Retrieve Force Density (N/m)
        const vector Fi_density = almF[segI];
        const vector Fip1_density = almF[segI + 1];

        // Normalize: (Density * Length) / Weight
        // This ensures the volume integral equals (Density * Length)
        const vector Fi_applied =
            Fi_density*ds/(nodeWeights[segI] + VSMALL);
        const vector Fip1_applied =
            Fip1_density*ds/(nodeWeights[segI + 1] + VSMALL);

        // Interpolation and Applying the force
        const vector Sj = -etaR*((1.0 - s)*Fi_applied + s*Fip1_applied);

        forceVals[c] = Sj;
        eqn.source()[c] += ((Sj/fluidRho_)*V[c]);
        ffvOption += ((Sj/fluidRho_)*V[c]);
    }
    reduce(ffvOption, sumOp<vector>());
    Info<< "Sum of ALM forces applied to Fluid (fvOptions) = "
        << ffvOption << endl;
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fvBeamPorosity::fvBeamPorosity
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    option(name, modelType, dict, mesh),
    forceFilePtr_(),
    eps_(),
    beamName_(),
    fluidRho_()
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
            postProcDir =
                mesh.time().path()
               /".."
               /"postProcessing"
               /startTimeName;
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

void fvBeamPorosity::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    applyBeamForce(eqn);
}


void fvBeamPorosity::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    applyBeamForce(eqn);
}


bool fvBeamPorosity::read(const dictionary& dict)
{
    if (option::read(dict))
    {
        if (!coeffs_.readIfPresent("UNames", fieldNames_))
        {
            fieldNames_.resize(1);
            fieldNames_.first() = coeffs_.getOrDefault<word>("U", "U");
        }
        option::resetApplied();
        coeffs_.readEntry("epsilon", eps_);
        coeffs_.readEntry("beamName", beamName_);
        coeffs_.readEntry("rho", fluidRho_);
        return true;
    }
    return false;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
