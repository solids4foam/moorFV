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

#include "sixDoFRigidBodyMotionRestraintNew.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::sixDoFRigidBodyMotionRestraintNew>
Foam::sixDoFRigidBodyMotionRestraintNew::New
(
    const dictionary& sDoFRBMRDict, const Time& time
)
{
    word sixDoFRigidBodyMotionRestraintNewTypeName =
        sDoFRBMRDict.lookup("sixDoFRigidBodyMotionRestraintNew");

    // Info<< "Selecting sixDoFRigidBodyMotionRestraintNew function "
    //     << sixDoFRigidBodyMotionRestraintNewTypeName << endl;

    dictionaryConstructorTable::iterator cstrIter =
    dictionaryConstructorTablePtr_->find
    (
        sixDoFRigidBodyMotionRestraintNewTypeName
    );

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "sixDoFRigidBodyMotionRestraintNew::New"
            "("
                "const dictionary& sDoFRBMRDict"
            ")"
        )   << "Unknown sixDoFRigidBodyMotionRestraintNew type "
            << sixDoFRigidBodyMotionRestraintNewTypeName << endl << endl
            << "Valid  sixDoFRigidBodyMotionRestraintNews are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<sixDoFRigidBodyMotionRestraintNew>
    (
        cstrIter()(sDoFRBMRDict, time)
    );
}


// ************************************************************************* //
