/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "finiteVolumeBeam.H"
#include "addToRunTimeSelectionTable.H"
#include "sixDOFODE.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(finiteVolumeBeam, 0);
    addToRunTimeSelectionTable
    (
        combinedRestraint,
        finiteVolumeBeam,
        word
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::finiteVolumeBeam::finiteVolumeBeam
(
    const word& name,
    const dictionary& dict,
    const sixDOFODE& sixDOF
)
:
    combinedRestraint(name, dict, sixDOF),
    beamPtr_()
{
    // Create beam model
    // Why does the beamModel need non-const access to time?
    beamPtr_.set
    (
        new beamModels::coupledTotalLagNewtonRaphsonBeam
        (
            const_cast<Time&>(sixDOF.dict().time())
        )
    );

    // Do we need to do anything else here to initialise the beam correctly?
}


Foam::autoPtr<Foam::combinedRestraint>
Foam::finiteVolumeBeam::clone() const
{
    return autoPtr<combinedRestraint>
    (
        new finiteVolumeBeam(*this)
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::finiteVolumeBeam::~finiteVolumeBeam()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::vector Foam::finiteVolumeBeam::restrainingForce
(
    const scalar t,
    const tensor& toRelative,
    const vector& x,
    const vector& u
) const
{
    // 1. Set the boundary conditions on the beam using the inputs above
    //    1a. Set the displacement condition at the attachment point using the
    //        new position 'x'
    //    1b. Set the rotation condition at the attachment point: we have two
    //        options (can be set at start; no need to update):
    //            - Option 1: set zero rotation
    //            - Option 2: set zero moment (seems more physical)

    // Lookup the displacement increment field
    volVectorField& DW = const_cast<volVectorField&>
    (
        sixDOF().dict().time().lookupObject<volVectorField>("DW")
    );

    // Find the attachment patch index
    const fvMesh& mesh = DW.mesh();
    const label patchID =
        mesh.boundaryMesh().findPatchID("name_of_beam_attachment_patch");

    if (patchID == -1)
    {
        FatalErrorIn("Foam::finiteVolumeBeam::restrainingForce(...)")
            << "Attachment patch not found!" << abort(FatalError);
    }

    // Set the value of DW on the patch
    // To-do: calculate what this should be as a function of the inputs above
    // and maybe mesh.boundary()[patchID].Cf()[0] and maybe use W.oldTime()
    DW.boundaryField()[patchID] == vector(1, 2, 3);

    // 2. Call the beam evolve function to solve the beam equations, i.e.
    beamPtr_().evolve();

    // 3. Retrieve and return the force at the attachment point
    // Lookup some force or stress field available within the beam model and
    // extract the value on the attachment patch, e.g. use lookupObject for the
    // force (maybe it is the surfaceVectorField "Q")
    const vector force = vector(1, 2, 3);

    WarningIn("restrainingForce")
        << "Force set to zero" << endl;

    return force;
}


Foam::vector Foam::finiteVolumeBeam::restrainingMoment
(
    const scalar t,
    const tensor& toRelative,
    const vector& omega
) const
{
    WarningIn("restrainingMoment")
        << "Moment set to zero" << endl;

    return vector::zero;
}


void Foam::finiteVolumeBeam::write(Ostream& os) const
{
    os.writeKeyword("type") << tab << type()
        << token::END_STATEMENT << nl << nl;

    // os.writeKeyword("coeff") << tab << coeff_
    //    << token::END_STATEMENT << nl;
}


// ************************************************************************* //
