/*---------------------------------------------------------------------------*\
License
    This file is part of moorFV.

    moorFV is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    moorFV is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with moorFV.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "sixDoFRigidBodyBeamMotionRestraint.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::sixDoFRigidBodyBeamMotionRestraint>
Foam::sixDoFRigidBodyBeamMotionRestraint::New
(
    const word& name,
    const dictionary& dict,
    const Time& time
)
{
    const word modelType
    (
        dict.get<word>("sixDoFRigidBodyBeamMotionRestraint")
    );

    auto cstrIter = dictionaryConstructorTablePtr_->cfind(modelType);

    if (!cstrIter.found())
    {
        FatalIOErrorInLookup
        (
            dict,
            "sixDoFRigidBodyBeamMotionRestraint",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<sixDoFRigidBodyBeamMotionRestraint>
    (
        cstrIter()(name, dict, time)
    );
}


// ************************************************************************* //
