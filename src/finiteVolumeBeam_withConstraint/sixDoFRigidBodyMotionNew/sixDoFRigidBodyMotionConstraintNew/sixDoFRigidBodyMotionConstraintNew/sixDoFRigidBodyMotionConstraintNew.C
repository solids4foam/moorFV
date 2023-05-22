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

#include "sixDoFRigidBodyMotionConstraintNew.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::sixDoFRigidBodyMotionConstraintNew, 0);

defineRunTimeSelectionTable(Foam::sixDoFRigidBodyMotionConstraintNew, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionConstraintNew::sixDoFRigidBodyMotionConstraintNew
(
    const dictionary& sDoFRBMCDict, const Time& time
)
:
    sDoFRBMCCoeffs_
    (
        sDoFRBMCDict.subDict
        (
            word(sDoFRBMCDict.lookup("sixDoFRigidBodyMotionConstraintNew"))
          + "Coeffs"
        )
    ),
    tolerance_(readScalar(sDoFRBMCDict.lookup("tolerance"))),
    relaxationFactor_
    (
        sDoFRBMCDict.lookupOrDefault<scalar>("relaxationFactor", 1)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionConstraintNew::~sixDoFRigidBodyMotionConstraintNew()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::sixDoFRigidBodyMotionConstraintNew::read
(
    const dictionary& sDoFRBMCDict
)
{
    tolerance_ = (readScalar(sDoFRBMCDict.lookup("tolerance")));

    relaxationFactor_ = sDoFRBMCDict.lookupOrDefault<scalar>
    (
        "relaxationFactor",
        1
    );

    sDoFRBMCCoeffs_ = sDoFRBMCDict.subDict(type() + "Coeffs");

    return true;
}


void Foam::sixDoFRigidBodyMotionConstraintNew::write(Ostream& os) const
{
    os.writeKeyword("tolerance")
        << tolerance_ << token::END_STATEMENT << nl;

    os.writeKeyword("relaxationFactor")
        << relaxationFactor_ << token::END_STATEMENT << nl;
}

// ************************************************************************* //
