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
    combinedRestraint(name, dict, sixDOF) //,
    // coeff_(dict.lookup("coeff"))
{}


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
    const scalar,
    const tensor&,
    const vector& x,
    const vector& u
) const
{
    //return - (linSpringCoeffs_ & x) - (linDampingCoeffs_ & u);
    WarningIn("restrainingMoment")
        << "Force set to zero" << endl;

    return vector::zero;
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
